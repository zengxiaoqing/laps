cdis    Forecast Systems Laboratory
cdis    NOAA/OAR/ERL/FSL
cdis    325 Broadway cdis    Boulder, CO     80303 
cdis 
cdis    Forecast Research Division cdis    Local Analysis and Prediction Branch 
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

        subroutine lapswind_plot(c_display,i4time_ref,lun,NX_L,NY_L,
     1                           NZ_L, MAX_RADARS,r_missing_data,
     1                           laps_cycle_time,maxstns)

!       1995        Steve Albers         Original Version
!       1995 Dec 8  Steve Albers         Automated pressure range
!       97-Aug-14     Ken Dritz     Added NX_L, NY_L, NZ_L as dummy args
!       97-Aug-14     Ken Dritz     Added MAX_RADARS as dummy arg
!       97-Aug-14     Ken Dritz     Added r_missing_data as dummy arg
!       97-Aug-14     Ken Dritz     Added laps_cycle_time as dummy arg
!       97-Aug-14     Ken Dritz     Added maxstns as dummy arg
!       97-Aug-14     Ken Dritz     Changed LAPS_DOMAIN_FILE to 'nest7grid'
!       97-Aug-14     Ken Dritz     Removed include of lapsparms.for
!       97-Aug-14     Ken Dritz     Pass NX_L, NY_L, r_missing_data, and
!                                   laps_cycle_time to plot_cont
!       97-Aug-14     Ken Dritz     Pass NX_L, NY_L (a second time) and
!                                   r_missing_data, laps_cycle_time to
!                                   plot_barbs
!       97-Aug-14     Ken Dritz     Pass NX_L, NY_L, and laps_cycle_time to
!                                   plot_grid
!       97-Aug-14     Ken Dritz     Pass NX_L, NY_L, and laps_cycle_time to
!                                   plot_cldpcp_type
!       97-Aug-14     Ken Dritz     Pass NX_L, NY_L, and laps_cycle_time to
!                                   plot_stations
!       97-Aug-14     Ken Dritz     Pass maxstns to plot_station_locations
!       97-Aug-17     Ken Dritz     Pass r_missing_data to divergence
!       97-Aug-25     Steve Albers  Removed equivalence for uv_2d.
!                                   Removed equivalence for slwc_int.
!                                   Removed equivalence for slwc_2d.
!                                   Removed /lapsplot_cmn1/ and /lapsplot_cmn2/
!       97-Sep-24     John Smart    Added display funtionality for
!                                   polar orbiter (lrs).
!       98-Mar-23        "          Added lvd subdirectory flexibility.

        real*4 lat(NX_L,NY_L),lon(NX_L,NY_L),topo(NX_L,NY_L)
        real*4 rlaps_land_frac(NX_L,NY_L)

        character*1 c_display
        character*1 cansw
        character*13 filename,a13_time
        character*3 c3_site
        character*4 c4_string
        character*5 c5_string
        character*4 c4_log
        character*7 c7_string
        character*9 c9_string,a9_start,a9_end
        character infile*70

        character i4_to_byte

        real clow,chigh,cint_ref
        data clow/-200./,chigh/+400/,cint_ref/10./

        integer*4 idum1_array(NX_L,NY_L)

        real*4 dum1_array(NX_L,NY_L)
        real*4 dum2_array(NX_L,NY_L)
        real*4 dum3_array(NX_L,NY_L)
        real*4 dum4_array(NX_L,NY_L)

      ! Used for "Potential" Precip Type
!       logical l_mask(NX_L,NY_L)
        logical iflag_mvd,iflag_icing_index,iflag_cloud_type
     1         ,iflag_bogus_w
        logical iflag_snow_potential

        integer*4 ibase_array(NX_L,NY_L)
        integer*4 itop_array(NX_L,NY_L)

        character*2 c_field,c_metacode,c_type
        character*33 c33_label

!       integer*4 ity,ily,istatus
!       data ity/35/,ily/1010/

        real*4 mspkt
        data mspkt/.518/

!       Stuff to read in WIND file
        integer*4 KWND
        parameter (KWND = 3)
        real*4 u_2d(NX_L,NY_L) ! WRT True North
        real*4 v_2d(NX_L,NY_L) ! WRT True North
        real*4 w_2d(NX_L,NY_L)
        real*4 liw(NX_L,NY_L)
        real*4 helicity(NX_L,NY_L)
        real*4 vas(NX_L,NY_L)
        real*4 cint
        real*4 uv_2d(NX_L,NY_L,2)
!       equivalence (uv_2d(1,1,1),u_2d),(uv_2d(1,1,2),v_2d)

        real*4 div(NX_L,NY_L)
        real*4 dir(NX_L,NY_L)
        real*4 spds(NX_L,NY_L)
        real*4 umean(NX_L,NY_L) ! WRT True North
        real*4 vmean(NX_L,NY_L) ! WRT True North

        real*4 sndr_po(19,NX_L,NY_L)

        character*3 var_2d
        character*150  directory
        character*31  ext
        character*10  units_2d
        character*125 comment_2d
        character*9 comment_a,comment_b

!       For reading in radar data
        real*4 dummy_array(NX_L,NY_L)
        real*4 radar_array(NX_L,NY_L)
        real*4 radar_array_adv(NX_L,NY_L)

        real*4 v_nyquist_in_a(MAX_RADARS)
        real*4 rlat_radar_a(MAX_RADARS), rlon_radar_a(MAX_RADARS) 
        real*4 rheight_radar_a(MAX_RADARS)
        integer*4 i4time_radar_a(MAX_RADARS)
        integer*4 n_vel_grids_a(MAX_RADARS)
        character*4 radar_name,radar_name_a(MAX_RADARS)
        character*31 ext_radar_a(MAX_RADARS)

        real*4 u_3d(NX_L,NY_L,NZ_L) ! WRT True North
        real*4 v_3d(NX_L,NY_L,NZ_L) ! WRT True North
        real*4 omega_3d(NX_L,NY_L,NZ_L)
        real*4 grid_ra_ref(NX_L,NY_L,NZ_L)
        real*4 grid_ra_vel(NX_L,NY_L,NZ_L,MAX_RADARS)
        real*4 grid_ra_nyq(NX_L,NY_L,NZ_L,MAX_RADARS)
!       real*4 grid_ra_rfill(NX_L,NY_L,NZ_L)

        real*4 lifted(NX_L,NY_L)
        real*4 height_2d(NX_L,NY_L)
        real*4 temp_2d(NX_L,NY_L)
        real*4 tw_sfc_k(NX_L,NY_L)
        real*4 td_2d(NX_L,NY_L)
        real*4 pres_2d(NX_L,NY_L)
        real*4 temp_3d(NX_L,NY_L,NZ_L)
        real*4 temp_col_max(NX_L,NY_L)
        real*4 rh_3d(NX_L,NY_L,NZ_L)
        real*4 pressures_mb(NZ_L)
        real*4 q_3d(NX_L,NY_L,NZ_L)
        real*4 slwc_3d(NX_L,NY_L,NZ_L)
        real*4 cice_3d(NX_L,NY_L,NZ_L)

!       real*4 slwc_int(NX_L,NY_L)
        real*4 column_max(NX_L,NY_L)
        character pcp_type_2d(NX_L,NY_L)
        character b_array(NX_L,NY_L)

        real*4 slwc_2d(NX_L,NY_L)
        real*4 cice_2d(NX_L,NY_L)
        real*4 field_2d(NX_L,NY_L)

        real*4 snow_2d(NX_L,NY_L)
        real*4 snow_2d_buf(NX_L,NY_L)
        real*4 precip_2d(NX_L,NY_L)
        real*4 precip_2d_buf(NX_L,NY_L)
        real*4 accum_2d(NX_L,NY_L)
        real*4 accum_2d_buf(NX_L,NY_L)

!       Local variables used in
        real*4 snow_accum_pd(NX_L,NY_L)
        real*4 snow_rate(NX_L,NY_L) ! M/S
        real*4 precip_rate(NX_L,NY_L) ! M/S
        real*4 dbz_2d(NX_L,NY_L)
        real*4 temp_sfc_k(NX_L,NY_L)
        real*4 td_sfc_k(NX_L,NY_L)
        real*4 pres_sta_pa(NX_L,NY_L)
        logical l_mask(NX_L,NY_L)
        integer ipcp_1d(NZ_L)


        character cldpcp_type_3d(NX_L,NY_L,NZ_L)
        real*4 mvd_3d(NX_L,NY_L,NZ_L)
        integer*4 iarg

        real*4 cloud_cvr(NX_L,NY_L)
        real*4 cloud_ceil(NX_L,NY_L)
        real*4 cloud_low(NX_L,NY_L)
        real*4 cloud_top(NX_L,NY_L)

        character asc9_tim_3dw*9, asc24_tim_3dw*24
        character asc9_tim_r*9, asc9_tim*9, asc_tim_9
        character asc9_tim_t*9
        character asc9_tim_n*9
        character c9_radarage*9

        character*255 c_filespec_ra
        character*255 c_filespec_src
        data c_filespec_ra /'/data/laps/nest7grid/lapsprd/vrc/*.vrc'/       
        data c_filespec_src/'*.src'/

        character*255 c_filespec
        character*255 cfname
        character*2   cchan

        logical lapsplot_pregen,l_precip_pregen,l_pregen,l_radar_read
        data lapsplot_pregen /.true./

        real*4 heights_3d(NX_L,NY_L,NZ_L)

        real*4 p_1d_pa(NZ_L)
        real*4 pbe_2d(NX_L,NY_L)
        real*4 nbe_2d(NX_L,NY_L)

        include 'laps_cloud.inc'

        real*4 clouds_3d(NX_L,NY_L,KCLOUD)
        real*4 cld_pres(KCLOUD)

        character*80 grid_fnam_common
        common / grid_fnam_cmn / grid_fnam_common

!       common/lapsplot_cmn1/heights_3d,slwc_3d,cice_3d,q_3d,rh_3d,temp_
!    13d,
!    1              u_3d,v_3d,omega_3d,grid_ra_ref,grid_ra_vel

        common /supmp1/ dummy,part

!       COMMON /CONRE1/IOFFP,SPVAL,EPSVAL,CNTMIN,CNTMAX,CNTINT,IOFFM

        include 'satellite_dims_lvd.inc'
        include 'satellite_common_lvd.inc'

        data mode_lwc/2/

        i_overlay = 0
        jdot = 1   ! 1 = Dotted County Boundaries, 0 = Solid County Boundaries
        part = 0.9 ! For plotting routines

        ioffm = 1 ! Don't plot label stuff in conrec

        grid_fnam_common = 'nest7grid'  ! Used in get_directory to modify
                                      ! extension based on the grid domain
        ext = 'nest7grid'

!       Get the location of the static grid directory
        call get_directory(ext,directory,len_dir)

        var_2d='LAT'
        call rd_laps_static (directory,ext,nx_l,ny_l,1,var_2d,
     1units_2d,comment_2d,
     1lat,rspacing_dum,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error reading LAPS static-lat'
            return
        endif

        var_2d='LON'
        call rd_laps_static (directory,ext,nx_l,ny_l,1,var_2d,
     1units_2d,comment_2d,
     1lon,rspacing_dum,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error reading LAPS static-lon'
            return
        endif

        var_2d='AVG'
        call rd_laps_static (directory,ext,nx_l,ny_l,1,var_2d,
     1units_2d,comment_2d,
     1topo,rspacing_dum,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error reading LAPS static-topo'
            return
        endif

        var_2d='LDF'
        call rd_laps_static (directory,ext,nx_l,ny_l,1,var_2d,
     1units_2d,comment_2d,
     1rlaps_land_frac,rspacing_dum,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error reading LAPS static-ldf'
            return
        endif

        if(lun .eq. 5)call logit('nest7grid')

1200    write(6,11)
11      format(//'  SELECT FIELD:  '
     1       /'     [wd] Wind, [wb,wr] Background Wind, [co] Cloud Omega
     1'
     1       /'     [lw] li*w, [li] li, [he] helicity, [pe] CAPE, [ne] C
     1AP'
     1       /'     [ra] Radar Data - NOWRAD vrc files,  [rx] Max Radar'
     1       /'     [rd] Radar Data - Doppler Ref-Vel (v01-v02...)'
     1       /
     1       /'     SFC: [p,pm,ps,tt,td,ws,vv,hu,ta,th,te,vo,mr,mc,dv'       
     1       ,',ha,ma,sp,cs,vs,tw,fw,hi]'
     1       /'          [ob,st] obs plot/station locations'
     1       /
     1       /'     [t]  Temperature, [pt] Potential Temperature,'
     1       ,'  [hh] Height of Const Temp Sfc'
     1       /
     1       /4x,' [ci] Cloud Ice     ',22x,'  [ls] Smith - Feddes Cloud
     1 LWC'
     1       /'     [is] Smith - Feddes Integrated Cloud LWC  '
     1       /'     [mv] Mean Volume Drop Diam,   [ic] Icing Index,'
     1       ,' [tc,tp] Cloud/Precip Type'
     1       /
     1       /'     [cc] Cld Ceiling (AGL), [cb] Lowest Cld Base (MSL)'
     1       ,',  [ct] Highest Cld Top'
     1       /'     [cv] Cloud Cover (2-D), [pw] Precipitable Water'
     1       /
     1       /'     [ht,hb,hy] Hts (LAPS,Bkgnd,Hydrstc),'
     1       ,' [sa/pa] Snow/Pcp Accum'
     1       /'     [sc] Snow Cover'
     1       /'     [sh,rh] Specific/Rel Humidity'
     1       ,'     [tr,lf,gr] Terrain/Land Frac/Grid, '
     1       /'     [v1,v2,v3,v4,v5,po] IR Twm/av; VCF/VIS; Tsfc-11u'
     1        '; Polar Orbiter'
     1       //' ',52x,'[q] quit/display ? '$)
        read(lun,12)c_type
12      format(a2)

        c4_log = 'h '//c_type
        if(lun .eq. 5 .and. c_type .ne. 'q ')call logit(c4_log)

!  ***  Ask for wind field ! ***************************************************
        if(    c_type .eq. 'wd' .or. c_type .eq. 'wb'
     1    .or. c_type .eq. 'co' .or. c_type .eq. 'wr')then

            if(c_type .eq. 'wd')then
                ext = 'lw3'
                call get_directory(ext,directory,len_dir)
                c_filespec = directory(1:len_dir)//'*.'//ext(1:3)
            elseif(c_type .eq. 'wb')then
                call make_fnam_lp(i4time_ref,asc9_tim_3dw,istatus)

                ext = 'lga'

                call get_directory(ext,directory,len_dir)
                c_filespec = directory(1:len_dir)//'*.'//ext(1:3)

            elseif(c_type .eq. 'wr')then
                call make_fnam_lp(i4time_ref,asc9_tim_3dw,istatus)

                ext = 'ram'

                call get_directory(ext,directory,len_dir)
                c_filespec = directory(1:len_dir)//'*.'//ext(1:3)

            elseif(c_type .eq. 'co')then
                ext = 'lco'
                call get_directory(ext,directory,len_dir)
                c_filespec = directory(1:len_dir)//'*.'//ext(1:3)
            endif

            write(6,*)
            write(6,*)'    Looking for 3D laps wind data:'
            call get_file_time(c_filespec,i4time_ref,i4time_3dw)
            call make_fnam_lp(I4time_3dw,asc9_tim_3dw,istatus)

            if(c_type .eq. 'wd')then
                write(6,13)
13              format(
     1    '     Enter Level in mb, -1 = steering, 0 = sfc',24x,'? '$)
            else
                write(6,14)
14              format('     Enter Level in mb',48x,'? '$)
            endif

            read(lun,*)k_level

            k_mb = k_level

            if(k_level .gt. 0)then
               k_level = nint(zcoord_of_pressure(float(k_level*100)))
            endif

            if(c_type .eq. 'co')then
                c_field = 'w'
                goto115
            endif

            if(k_level .gt. -1)then

                if(k_level .eq. 0)then ! SFC Winds
                    write(6,102)
102                 format(/
     1 '  Field [di,sp,u,v,vc (barbs), ob (obs)]',30x,'? '$)
                    read(lun,12)c_field
                    ext = 'lwm'
                    call get_directory(ext,directory,len_dir)
                    c_filespec = directory(1:len_dir)//'*.'//ext(1:3)

                    var_2d = 'SU'
                    call get_laps_2d(i4time_3dw,ext,var_2d
     1          ,units_2d,comment_2d,NX_L,NY_L,u_2d,istatus)
                    var_2d = 'SV'
                    call get_laps_2d(i4time_3dw,ext,var_2d
     1          ,units_2d,comment_2d,NX_L,NY_L,v_2d,istatus)


                else if(k_level .gt. 0)then
                    write(6,103)
103                 format(/
     1 '  Field [di,sp,u,v,w,dv,vc (barbs), ob (obs))]'
     1                                          ,24x,'? '$)
                    read(lun,12)c_field
                    write(6,*)' ext = ',ext
                    if(c_field .ne. 'w ')then
                      call get_uv_2d(i4time_3dw,k_level,uv_2d,
     1                                          ext,NX_L,NY_L,istatus)
                      call move(uv_2d(1,1,1),u_2d,NX_L,NY_L)
                      call move(uv_2d(1,1,2),v_2d,NX_L,NY_L)
                    endif

                endif


            elseif(k_level .eq. -1)then ! Calculate mean winds from 3d grids

                write(6,*)' Getting pregenerated mean wind file'

                ext = 'lwm'
                call get_directory(ext,directory,len_dir)
                var_2d = 'MU'
                call get_laps_2d(i4time_3dw,ext,var_2d
     1          ,units_2d,comment_2d,NX_L,NY_L,u_2d,istatus)
                var_2d = 'MV'
                call get_laps_2d(i4time_3dw,ext,var_2d
     1          ,units_2d,comment_2d,NX_L,NY_L,v_2d,istatus)

                write(6,104)
104             format(/
     1     '  Field [di,sp,u,v,vc (barbs)]   ',25x,'? '$)
                read(lun,12)c_field
            endif

!  ***      Display Wind Data  ******************************************************

115         if(c_field .eq. 'di' .or. c_field .eq. 'sp')then
                do i = 1,NX_L
                do j = 1,NY_L
                    if(u_2d(i,j) .eq. r_missing_data
     1    .or. v_2d(i,j) .eq. r_missing_data)then
                        dir(i,j)  = r_missing_data
                        spds(i,j) = r_missing_data
                    else
                        call uvgrid_to_disptrue(u_2d(i,j),
     1                                  v_2d(i,j),
     1                                  dir(i,j),
     1                                  spds(i,j),
     1                                  lon(i,j)     )
                        spds(i,j) = spds(i,j) / mspkt
                    endif
                enddo ! j
                enddo ! i
            endif

            if(c_field .eq. 'di')then
                call mklabel33
     1      (k_level,' Isogons   (deg)   ',c33_label)

                call plot_cont(dir,1e0,clow,chigh,30.,asc9_tim_3dw,
     1  c33_label,i_overlay,c_display,'nest7grid',lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

            else if(c_field .eq. 'sp')then
                call mklabel33
     1              (k_level,' Isotachs (kt)     ',c33_label)
                if(k_level .gt. 0 .and. k_mb .le. 500)then
                    cint = 10.
                else
                    cint = 5.
                endif

                call plot_cont(spds,1e0,0.,300.,cint,asc9_tim_3dw,
     1  c33_label,i_overlay,c_display,'nest7grid',lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

            else if(c_field .eq. 'u ')then
                call mklabel33
     1      (k_level,' U - Component      ',c33_label)

                call plot_cont(u_2d,1e0,clow,chigh,10.,asc9_tim_3dw,
     1  c33_label,i_overlay,c_display,'nest7grid',lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

            else if(c_field .eq. 'v ')then
                call mklabel33
     1      (k_level,' V - Component      ',c33_label)

                call plot_cont(v_2d,1e0,clow,chigh,10.,asc9_tim_3dw,
     1  c33_label,i_overlay,c_display,'nest7grid',lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

            else if(c_field .eq. 'vc' .or. c_field .eq. 'ob')then
                call mklabel33
     1      (k_level,' WINDS    (kt)     ',c33_label)

                if(max(NX_L,NY_L) .gt. 50)then
                    interval = 2
                else
                    interval = 1
                endif

                size = float(interval) * .14

                call plot_barbs(u_2d,v_2d,lat,lon,topo,size,interval
     1               ,asc9_tim_3dw
     1               ,c33_label,c_field,k_level,i_overlay,c_display       
     1               ,'nest7grid',NX_L,NY_L,NZ_L,grid_ra_ref,grid_ra_vel      
     1               ,NX_L,NY_L,r_missing_data,laps_cycle_time,jdot)

            else if(c_field .eq. 'w')then ! Display W fields
!               if(lapsplot_pregen .and. k_level .eq. 7)then
                if(.false.)then
                    write(6,*)' Getting pregenerated 600mb omega file'

                    var_2d = 'OM'
                    ext = 'liw'
                    call get_laps_2d(i4time_3dw,ext,var_2d
     1                  ,units_2d,comment_2d,NX_L,NY_L,w_2d,istatus)

                else
                    var_2d = 'OM'
                    ext = 'lw3'
                    call get_laps_2dgrid(i4time_3dw,0,i4time_nearest,
     1              ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                          ,w_2d,k_mb,istatus)

                endif

                call mklabel33
     1      (k_level,' Omega    (ubar/s) ',c33_label)

                do j = 1,NY_L
                do i = 1,NX_L
                    if(w_2d(i,j) .eq. r_missing_data)then
                        w_2d(i,j) = 0.
                    endif
                enddo ! i
                enddo ! j

                call plot_cont(w_2d,1e-1,0.,0.,-1.0,asc9_tim_3dw,
     1  c33_label,i_overlay,c_display,'nest7grid',lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

            elseif(c_field .eq. 'dv')then ! Display Divergence Field
                call divergence(u_2d,v_2d,div,lat,lon,NX_L,NY_L
     1  ,dum1_array,dum2_array
     1  ,dum3_array,dum4_array,dummy_array,radar_array,r_missing_data)
                call mklabel33(k_level,' DVRGNC  1e-5 s(-1)',c33_label)

                call plot_cont(div,1e-5,0.,0.,0.,asc9_tim_3dw,
     1  c33_label,i_overlay,c_display,'nest7grid',lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

            endif

        elseif(c_type .eq. 'lw')then ! Read in Liw field from 3d grids

            write(6,*)
            write(6,*)'    Looking for laps li*w data:'

            if(lapsplot_pregen)then
                ext = 'liw'
                call get_directory(ext,directory,len_dir)
                c_filespec = directory(1:len_dir)//'*.'//ext(1:3)
                call get_file_time(c_filespec,i4time_ref,i4time_3dw)
                call make_fnam_lp(I4time_3dw,asc9_tim_3dw,istatus)

                write(6,*)' Getting pregenerated Li * omega file'

                var_2d = 'LIW'
                ext = 'liw'
                call get_laps_2d(i4time_3dw,ext,var_2d
     1          ,units_2d,comment_2d,NX_L,NY_L,liw,istatus) ! K-Pa/s

            else ! Calculate LI * omega on the fly
                write(6,*)'    Looking for 3D laps wind data:'
                ext = 'lw3'
                call get_directory(ext,directory,len_dir)
                c_filespec = directory(1:len_dir)//'*.'//ext(1:3)
                call get_file_time(c_filespec,i4time_ref,i4time_3dw)
                call make_fnam_lp(I4time_3dw,asc9_tim_3dw,istatus)

                var_2d = 'OM'
                ext = 'lw3'
                lvl_2d = 600
                call get_laps_2dgrid(i4time_3dw,0,i4time_nearest,
     1              ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                          ,w_2d,lvl_2d,istatus)

!               Read in LI data
                var_2d = 'LI'
                ext = 'lsx'
                call get_laps_2dgrid(i4time_3dw,laps_cycle_time,i4time_n
     1earest,
     1          ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                          ,lifted,0,istatus)

                if(istatus .ne. 1)then
                    write(6,*)' Error reading Lifted Index data'
                    stop
                endif

                call cpt_liw(lifted,w_2d,NX_L,NY_L,liw) ! K-Pa/s

            endif ! Pregenerated LI * omega field

!           Logarithmically scale the values for display
!           do j = 1,NY_L,1
!           do i = 1,NX_L,1

!               if(liw(i,j) .ge. 3.16)then
!                   liw(i,j) = alog10(liw(i,j))
!               elseif(liw(i,j) .ge. 0.)then
!                   liw(i,j) = liw(i,j) * 0.5/3.16
!               endif

!           enddo ! j
!           enddo ! i

            chigh = 50.

            call plot_cont(liw,1e0,0.0,50.0,-0.5,asc9_tim_3dw,
     1   'LAPS sfc LI X 600mb omega  Pa-K/s',i_overlay,c_display
     1          ,'nest7grid',lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

        elseif(c_type .eq. 'li')then ! Read in Li field from 3d grids
            if(lapsplot_pregen)then
                write(6,*)' No pregenerated li file present, getting li
     1fm LSX'
!               Read in LI data
                var_2d = 'LI'
                ext = 'lsx'
                call get_laps_2dgrid(i4time_ref,7200,i4time_nearest,
     1          ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                          ,lifted,0,istatus)

            else
                call get_laps_2dgrid(i4time_ref,7200,i4time_nearest,
     1          ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                          ,lifted,0,istatus)

            endif

            if(istatus .ne. 1)then
                write(6,*)' Error reading Lifted Index data'
                goto1200
            endif

            call make_fnam_lp(i4time_nearest,asc9_tim_n,istatus)

            call plot_cont(lifted,1e-0,-20.,+40.,2.,asc9_tim_n,
     1      'LAPS    SFC Lifted Index     (K) ',i_overlay
     1  ,c_display,'nest7grid',lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

        elseif(c_type .eq. 'tw')then
            i4time_temp = i4time_ref / laps_cycle_time * laps_cycle_time

!           Read in surface temp data
            var_2d = 'T'
            ext = 'lsx'
            call get_laps_2dgrid(i4time_ref,laps_cycle_time,i4time_temp,
     1      ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                          ,temp_2d,0,istatus)

            if(istatus .ne. 1)then
                write(6,*)' LAPS Sfc Temp not available'
                return
            endif

!           Read in surface dewpoint data
            var_2d = 'TD'
            ext = 'lsx'
            call get_laps_2d(i4time_ref,
     1  ext,var_2d,units_2d,comment_2d,
     1                          NX_L,NY_L,td_2d,istatus)

            if(istatus .ne. 1)then
                write(6,*)' LAPS Sfc Dewpoint not available'
                return
            endif

!           Read in surface pressure data
            var_2d = 'PS'
            ext = 'lsx'
            call get_laps_2d(i4time_ref,
     1  ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                          ,pres_2d,istatus)

            if(istatus .ne. 1)then
                write(6,*)' LAPS Sfc Pressure not available'
                return
            endif

            call get_tw_approx_2d(temp_2d,td_2d,pres_2d,NX_L,NY_L,tw_sfc
     1_k)

            call make_fnam_lp(i4time_temp,asc9_tim_n,istatus)

            zero_c = 273.15

            do i = 1,NX_L
            do j = 1,NY_L
                field_2d(i,j) = tw_sfc_k(i,j) - zero_c
            enddo ! j
            enddo ! i


            call plot_cont(field_2d,1e-0,-30.,+30.,2.,asc9_tim_n,
     1      'LAPS    SFC Wetbulb (approx) (C) '
     1  ,i_overlay,c_display,'nest7grid',lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

        elseif(c_type .eq. 'ms' .or. c_type .eq. 'ob'
     1                          .or. c_type .eq. 'st'   )then
            i4time_plot = i4time_ref ! / laps_cycle_time * laps_cycle_time
            call get_filespec('lso',2,c_filespec,istatus)
            call get_file_time(c_filespec,i4time_ref,i4time_plot)
            call make_fnam_lp(i4time_plot,asc_tim_9,istatus)

            if(c_type .eq. 'st')iflag = 0
            if(c_type .eq. 'ms')iflag = 1
            if(c_type .eq. 'ob')iflag = 2

            c33_label = '                                 '

            call plot_stations(asc_tim_9,c33_label,c_field,i_overlay
     1   ,c_display,lat,lon,c_file,iflag
     1   ,NX_L,NY_L,laps_cycle_time)

        elseif(c_type .eq. 'he')then
            write(6,*)
!           write(6,*)'    Looking for 3D laps wind data:'
!           call get_file_time(c_filespec,i4time_ref,i4time_3dw)

            if(lapsplot_pregen)then
                write(6,*)' Getting pregenerated helicity file'
                var_2d = 'LHE'
                ext = 'lhe'
!               call get_laps_2d(i4time_3dw,ext,var_2d
!       1               ,units_2d,comment_2d,NX_L,NY_L,helicity,istatus)

                call get_laps_2dgrid(i4time_ref,7200,i4time_nearest,
     1          ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                  ,helicity,0,istatus)

                i4time_3dw = i4time_nearest
                call make_fnam_lp(I4time_3dw,asc9_tim_3dw,istatus)

            else
                ext = 'lw3'
                call get_uv_3d(i4time_3dw,NX_L,NY_L,NZ_L,u_3d,v_3d,ext,i
     1status)
!               call mean_wind(u_3d,v_3d,topo,NX_L,NY_L,NZ_L
!       1               ,dum1_array,dum2_array,dum3_array,idum1_array
!       1                                                       ,umean,vmean
!       1                                               ,u_2d,v_2d,istatus)
                call helicity_laps(u_3d,v_3d,u_2d,v_2d,topo
     1          ,dum1_array,dum2_array,dum3_array,idum1_array
     1          ,NX_L,NY_L,NZ_L
     1                          ,helicity,istatus)

                i4time_3dw = i4time_ref
                call make_fnam_lp(I4time_3dw,asc9_tim_3dw,istatus)

            endif

            c_field = 'he'
            kwind = 0
            clow = -100.
            chigh = +100.
            call plot_cont(helicity,1e-4,clow,chigh,5.,asc9_tim_3dw,
     1    'LAPS Helicity  sfc-500mb 1e-4m/s**2'
     1  ,i_overlay,c_display,'nest7grid',lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

        elseif(c_type .eq. 'v1' .or. c_type .eq. 'v2'
     1    .or. c_type .eq. 'v3' .or. c_type .eq. 'v4'
     1    .or. c_type .eq. 'v5')then
            write(6,*)
            write(6,*)'    Looking for Laps LVD data:'

            ext = 'lvd'

            if(iflag_lvd_common.ne.1)then
               call config_satellite_lvd(istatus)
               if(istatus.ne.1)then
                  return
               endif
            endif

            do k=1,maxsat
             if(isats(k).eq.1)then
              write(6,114)c_sat_id(k)
114           format(5x,'plot the data for ',a6,45x,' [y/n]? ',$)
              read(lun,*)cansw
              if(cansw.eq.'y'.or.cansw.eq.'Y')then

               call get_directory(ext,directory,len_dir)
               directory=directory(1:len_dir)//c_sat_id(k)//'/'

               if(c_type .eq. 'v1')then

                var_2d = 'S8W'
                call get_2dgrid_dname(directory
     1          ,i4time_ref,10000,i4time_nearest
     1          ,ext,var_2d
     1          ,units_2d,comment_2d,NX_L,NY_L,vas,0,istatus)

                if(istatus .ne. 1)then
                    write(6,*)' Cant find S8W Analysis ',istatus
                    goto1200
                endif
                c33_label = 'LAPS IR Skin Temps Deg C  -'//c_sat_id(k)

               elseif(c_type .eq. 'v2')then

                var_2d = 'S8A'
                call get_2dgrid_dname(directory
     1          ,i4time_ref,10000,i4time_nearest
     1          ,ext,var_2d
     1          ,units_2d,comment_2d,NX_L,NY_L,vas,0,istatus)

                if(istatus .ne. 1)then
                    write(6,*)' Cant find VAS/S8A Analysis'
                    goto1200
                endif
                c33_label = 'LAPS Ave 11.2 u     Deg C -'//c_sat_id(k)

               elseif(c_type .eq. 'v3')then

                var_2d = 'ALB'
                call get_2dgrid_dname(directory
     1          ,i4time_ref,10000,i4time_nearest
     1          ,ext,var_2d
     1          ,units_2d,comment_2d,NX_L,NY_L,vas,0,istatus)

                if(istatus .ne. 1 .and. istatus .ne. -1)then
                    write(6,*)' Cant find ALB Analysis'
                    goto1200
                endif
                c33_label = 'LAPS VIS Cloud Fraction   -'//c_sat_id(k)

                do i = 1,NX_L
                do j = 1,NY_L
                    vas(i,j) = albedo_to_cloudfrac(vas(i,j))
                enddo 
                enddo

               elseif(c_type .eq. 'v4')then

                var_2d = 'SVS'
                call get_2dgrid_dname(directory
     1          ,i4time_ref,10000,i4time_nearest
     1          ,ext,var_2d
     1          ,units_2d,comment_2d,NX_L,NY_L,vas,0,istatus)
                if(istatus .ne. 1 .and. istatus .ne. -1)then
                    write(6,*)' Cant find VIS Analysis'
                    goto1200
                endif
                c33_label = 'LAPS VIS (Unnormalized)   -'//c_sat_id(k)

               elseif(c_type .eq. 'v5')then

                var_2d = 'S8A'
                call get_2dgrid_dname(directory
     1          ,i4time_ref,10000,i4time_nearest
     1          ,ext,var_2d
     1          ,units_2d,comment_2d,NX_L,NY_L,vas,0,istatus)
                if(istatus .ne. 1)then
                    write(6,*)' Cant find VAS/S8A Analysis'
                    goto1200
                endif
                c33_label = 'LAPS SFC T - Band 8  (K)  -'//c_sat_id(k)

!               Get sfc T to take the difference...
                ext = 'lsx'
                var_2d = 'T'
                call get_laps_2dgrid(i4time_nearest,0,i4time_nearest
     1                ,ext,var_2d
     1                ,units_2d,comment_2d,NX_L,NY_L,dum1_array,0,istatu
     1s)
                if(istatus .ne. 1)then
                    write(6,*)' Cant find VAS/S8A Analysis'
                    goto1200
                endif

                do i = 1,NX_L
                do j = 1,NY_L
                    vas(i,j) = dum1_array(i,j) - vas(i,j)
                enddo
                enddo

               endif

               if(c_type .eq. 'v1' .or. c_type .eq. 'v2')then
                do i = 1,NX_L
                do j = 1,NY_L
                 vas(i,j) = vas(i,j) - 273.15
                enddo
                enddo
                clow = -80.
                chigh = +40.
                cint = 10.
               elseif(c_type .eq. 'v3')then
                clow = -0.5
                chigh = +1.5
                cint = 0.2
               elseif(c_type .eq. 'v4')then
                clow = 0.0
                chigh = 256.
                cint = 05.
               elseif(c_type .eq. 'v5')then
                clow = -8.0
                chigh = 20.
                cint = 4.
               endif

               call make_fnam_lp(i4time_nearest,asc9_tim,istatus)

               call plot_cont(vas,1e0,clow,chigh,cint,asc9_tim,
     1    c33_label,i_overlay,c_display,'nest7grid',lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

              endif
             endif
            enddo

        elseif( c_type .eq. 'po' )then

          call make_fnam_lp(i4time_ref,asc9_tim,istatus)
          ext = 'lsr'
          call get_directory(ext,directory,len_dir)
          cfname=directory(1:len_dir)//asc9_tim//'_12.lsr'
          open(14,file=cfname,form='unformatted',status='old',err=18)
          goto 27
18        cfname=directory(1:len_dir)//asc9_tim//'_14.lsr'
          open(14,file=cfname,form='unformatted',status='old',err=19)

          n=index(cfname,' ')-1
          write(6,*)'Reading ',cfname(1:n)
          goto 28
27        write(6,*)'Reading ',cfname(1:n)
28        read(14)sndr_po
          write(6,*)'Enter the channel [1-19]'
          read(5,25)ichan
25        format(i2)
          do j=1,NY_L
          do i=1,NX_L
             if(sndr_po(ichan,i,j).ne.r_missing_data)then
                vas(i,j)=sndr_po(ichan,i,j)-273.15
             endif
          enddo
          enddo

          clow = -80.
          chigh = +40.
          cint = 5.
          write(cchan,111)ichan
 111      format(i2)
          if(cchan(1:1).eq.' ')cchan(1:1)='0'
          if(cchan(2:2).eq.' ')cchan(2:2)='0'

          c33_label = 'LAPS Polar Orbiter Channel '//cchan//' deg C'

          call plot_cont(vas,1e0,clow,chigh,cint,asc9_tim,
     1    c33_label,i_overlay,c_display,'nest7grid',lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

          goto 21
19        write(6,*)'Not able to open an lsr file ', asc9_tim
21        continue

        elseif( c_type .eq. 'ra' .or. c_type .eq. 'gc'
     1    .or.  c_type .eq. 'rr'
     1    .or.  c_type .eq. 'rd'                          )then

            if(c_type .eq. 'ra')mode = 1
            if(c_type .eq. 'gc')mode = 2

            i4time_tmp1 = (i4time_ref)/laps_cycle_time * laps_cycle_time
            i4time_tmp2 = (i4time_ref-2400)/laps_cycle_time * laps_cycle
     1_time

            if(c_type .eq. 'rr')then
                if(i4time_ref .ne. i4time_tmp1)then
                    i4time_get = i4time_tmp2
                else
                    i4time_get = i4time_ref
                endif
            else
                i4time_get = i4time_ref
            endif

2010        if(.not. l_radar_read)then

              if(c_type .ne. 'rd')then ! Read data from vrc files

!               Obtain height field
                ext = 'lt1'
                var_2d = 'HT'
                call get_laps_3dgrid(
     1                   i4time_get,10000000,i4time_ht,
     1                   NX_L,NY_L,NZ_L,ext,var_2d
     1                  ,units_2d,comment_2d,heights_3d,istatus)
                if(istatus .ne. 1)then
                    write(6,*)' Error locating height field'
                    return
                endif

                call get_radar_ref(i4time_get,100000,i4time_radar,mode
     1            ,.true.,NX_L,NY_L,NZ_L,lat,lon,topo,.true.,.true.
     1            ,heights_3d
     1            ,grid_ra_ref,n_ref
     1            ,rlat_radar,rlon_radar,rheight_radar,istat_2dref
     1            ,istat_3dref)

              else ! 'rd' option: read data from v01, v02, etc.

                write(6,*)' Reading velocity data from the radars'

                call get_multiradar_vel(
     1            i4time_get,100000000,i4time_radar_a
     1           ,max_radars,n_radars,ext_radar_a,r_missing_data
     1           ,.true.,NX_L,NY_L,NZ_L
     1           ,grid_ra_vel,grid_ra_nyq,v_nyquist_in_a
     1           ,n_vel_grids_a
     1           ,rlat_radar_a,rlon_radar_a,rheight_radar_a,radar_name_a       
     1           ,istat_radar_vel,istat_radar_nyq)

                if(istat_radar_vel .eq. 1)then
                  write(6,*)' Radar 3d vel data successfully read in'
     1                       ,(n_vel_grids_a(i),i=1,n_radars)
                else
                  write(6,*)' Radar 3d vel data NOT successfully read in
     1'
     1                       ,(n_vel_grids_a(i),i=1,n_radars)
                  return
                endif


!               Ask which radar number (extension)
                write(6,*)
                write(6,2026)
2026            format('         Enter Radar # (for reflectivity)  '
     1                                                     ,27x,'? '$)
                read(lun,*)i_radar

                write(6,*)' Reading reflectivity data from radar '
     1                                                       ,i_radar

!               Obtain radar time (call get_file_time - ala put_derived_wind)

!               Obtain height field
                ext = 'lt1'
                var_2d = 'HT'
                call get_laps_3dgrid(
     1                   i4time_radar_a(i_radar),1000000,i4time_ht,
     1                   NX_L,NY_L,NZ_L,ext,var_2d
     1                  ,units_2d,comment_2d,heights_3d,istatus)

                write(6,*)

                call read_radar_3dref(i4time_radar_a(i_radar),
!    1               0,i4_dum
     1               .true.,NX_L,NY_L,NZ_L,ext_radar_a(i_radar),
     1               lat,lon,topo,.true.,.true.,
     1               heights_3d,
     1               grid_ra_ref,
     1               rlat_radar,rlon_radar,rheight_radar,radar_name,
     1               n_ref_grids,istat_radar_2dref,istat_radar_3dref)

                     i4time_radar = i4time_radar_a(i_radar)

              endif

              call make_fnam_lp(i4time_radar,asc9_tim_r,istatus)
              l_radar_read = .true.
            endif

2015        write(6,2020)
2020        format(/'  [ve] Velocity Contours, '  
     1     ,' [vi] Velocity Image (no map)'
     1     /'  [rf] Reflectivity Data, '
     1           /'  [mr] Max Reflectivity, [vl] VIL, [mt] Max Tops,'
     1     /'  [lr] Low Lvl Reflectivity, '
     1     ,'[f1] 1 HR Fcst Max Reflectivity,'
     1     /' ',61x,' [q] Quit ? '$)
            read(lun,12)c_field

            if(  c_field .eq. 'rf' 
     1      .or. c_field .eq. 'vi' .or. c_field .eq. 've')then
                write(6,2021)
2021            format('         Enter Level in mb ',45x,'? '$)
                read(lun,*)k_level

                if(k_level .gt. 0)then
                    k_level = nint(zcoord_of_pressure(float(k_level*100)
     1))
                endif

            endif

            if(c_field .eq. 'mr')then ! Reflectivity data

!               if(lapsplot_pregen)then
                if(.false.)then
                    write(6,*)' Getting pregenerated radar data file'

                    var_2d = 'R00'
                    ext = 'lmr'
                    i4time_hour = (i4time_radar+laps_cycle_time/2)
     1                          /laps_cycle_time * laps_cycle_time
                    call make_fnam_lp(i4time_hour,asc9_tim_r,istatus)
                    call get_laps_2d(i4time_hour,ext,var_2d
     1          ,units_2d,comment_2d,NX_L,NY_L,radar_array,istatus)

                else
                    call get_max_ref(grid_ra_ref,NX_L,NY_L,NZ_L,radar_ar
     1ray)
                    call make_fnam_lp(i4time_radar,asc9_tim_r       
     1                                                      ,istatus)

                endif

!               Display R field

                call plot_cont(radar_array,1e0,0.,chigh,cint_ref,
     1        asc9_tim_r,'LAPS  Column Max Reflectivity    ',
     1          i_overlay,c_display,'nest7grid',lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

            elseif(c_field .eq. 'lr')then ! Low Lvl Reflectivity data
                i4time_hour = (i4time_radar+laps_cycle_time/2)
     1                          /laps_cycle_time * laps_cycle_time

                var_2d = 'LLR'
                ext = 'lmt'
                call get_laps_2dgrid(i4time_hour,laps_cycle_time*100
     1                              ,i4time_lr,ext,var_2d,units_2d
     1                              ,comment_2d,NX_L,NY_L
     1                              ,radar_array,0,istatus)

                call make_fnam_lp(i4time_lr,asc9_tim_r,istatus)
!               call get__low__ref(grid_ra_ref,topo,NX_L,NY_L,NZ_L,radar_array)

!               Display R field
                call plot_cont(radar_array,1e0,0.,chigh,cint_ref
     1                        ,asc9_tim_r
     1                        ,'LAPS Low LVL Reflectivity   (DBZ)'
     1                        ,i_overlay,c_display,'nest7grid',lat,lon
     1                        ,jdot,NX_L,NY_L,r_missing_data
     1                        ,laps_cycle_time)

            elseif(c_field .eq. 've')then
                call mklabel33(k_level,'  Radial Vel  (kt) ',c33_label)

                write(6,2031)
2031            format('         Enter Radar # (of ones available)  '
     1                                                    ,28x,'? '$)
                read(lun,*)i_radar

                call make_fnam_lp(i4time_radar_a(i_radar),asc9_tim_r
     1                                                      ,istatus)

                call plot_cont(grid_ra_vel(1,1,k_level,i_radar),.518
     1                        ,clow,chigh,5.,asc9_tim_r,c33_label
     1                        ,i_overlay,c_display,'nest7grid',lat,lon        
     1                        ,jdot,NX_L,NY_L,r_missing_data
     1                        ,laps_cycle_time)                          

            elseif(c_field .eq. 'vi')then
                call mklabel33(k_level,'  Radial Vel  (kt) ',c33_label)

                write(6,2031)
                read(lun,*)i_radar

                call plot_obs(k_level,.false.,asc9_tim,i_radar
     1          ,NX_L,NY_L,NZ_L,grid_ra_ref,grid_ra_vel(1,1,1,i_radar)
     1          ,lat,lon,topo,2)

            elseif(c_field .eq. 'rf')then
                call mklabel33(k_level,'   Reflectivity    ',c33_label)

                call plot_cont(grid_ra_ref(1,1,k_level),1e0,0.,chigh
     1                        ,cint_ref,asc9_tim_r,c33_label,i_overlay
     1                        ,c_display,'nest7grid',lat,lon,jdot,NX_L        
     1                        ,NY_L,r_missing_data,laps_cycle_time)

            elseif(c_field .eq. 'vl')then ! Do VIL

!               Initialize Radar Array
                do i = 1,NX_L
                do j = 1,NY_L
                    radar_array(i,j) = 0.
                enddo
                enddo

                do i = 1,NX_L
                do j = 1,NY_L
                do k = 1,NZ_L
                    radar_array(i,j) =
     1          max(radar_array(i,j),grid_ra_ref(i,j,k))

                enddo
                enddo
                enddo

                write(6,*)' Calculating VIL'
!               Call VIL(grid_ra_ref,radar_array)

                call plot_cont(radar_array,1e0,clow,chigh,10.
     1                        ,asc9_tim_r
     1                        ,'LAPS DUMMY VIL                   '
     1                        ,i_overlay,c_display,'nest7grid',lat,lon
     1                        ,jdot
     1                        ,NX_L,NY_L,r_missing_data,laps_cycle_time)       

            elseif(c_field .eq. 'mt')then ! Do Max Tops
                i4time_hour = (i4time_radar+laps_cycle_time/2)
     1                          /laps_cycle_time * laps_cycle_time

                var_2d = 'LMT'
                ext = 'lmt'
                call get_laps_2dgrid(i4time_hour,laps_cycle_time*100
     1                              ,i4time_lr,ext,var_2d,units_2d
     1                              ,comment_2d,NX_L,NY_L
     1                              ,radar_array,0,istatus)

!               call get_maxtops(grid_ra_ref,NX_L,NY_L,NZ_L,radar_array)

                highest_top_m = 0.

                do i = 1,NX_L
                do j = 1,NY_L
                    highest_top_m = max(radar_array(i,j),highest_top_m)
                enddo ! j
                enddo ! i

!               Generate Contour Range and Interval
                cont_high = int(highest_top_m/1000.)
                cont_low = int(max(cont_high/2.,4.))

                if(cont_high - cont_low .gt. 4)then
                    cint = 2.
                else
                    cint = 1.
                endif

!               Create a floor to the array for better contouring
                do i = 1,NX_L
                do j = 1,NY_L
                    rfloor = cont_low - 0.1
                    radar_array(i,j) = radar_array(i,j) / 1000.
                    radar_array(i,j) =
     1                  max(radar_array(i,j),rfloor)
                enddo ! j
                enddo ! i

!               Display Max Tops
                write(6,*)' Displaying Max Tops, cint = '
     1                          ,cont_low,cont_high,cint
                call plot_cont(radar_array,1e0,cont_low
     1                        ,cont_high,cint,asc9_tim_r
     1                        ,'LAPS Max Echo Tops    (km MSL)   '
     1                        ,i_overlay,c_display,'nest7grid',lat,lon
     1                        ,jdot,NX_L,NY_L,r_missing_data
     1                        ,laps_cycle_time)

            elseif(c_field .eq. 'f1')then ! Fcst Max Reflectivity

                if(lapsplot_pregen)then
                    write(6,*)' Getting pregenerated radar fcst file'

                    var_2d = 'R06'
                    ext = 'lmr'
                    i4time_hour = (i4time_radar+laps_cycle_time/2)
     1                          /laps_cycle_time * laps_cycle_time
                    call make_fnam_lp(i4time_hour,asc9_tim_r,istatus)
                    call get_laps_2d(i4time_hour,ext,var_2d
     1         ,units_2d,comment_2d,NX_L,NY_L,radar_array_adv,istatus)

                else

                    call get_max_ref(grid_ra_ref,NX_L,NY_L,NZ_L,radar_ar
     1ray)

                    kmax = nint(height_to_zcoord(5000.,istatus))


!                   Match the i4time of the wind analysis to the Radar Data
                    write(6,*)
                    write(6,*)'    Looking for 3D laps wind data:'
                    call get_file_time
     1          (c_filespec,i4time_radar,i4time_3dw)
                    call make_fnam_lp(I4time_3dw,asc9_tim_3dw,istatus)

                    write(6,*)' Using Latest 3D laps data at '
     1                           ,Asc9_tim_3dw
                    call get_uv_3d(i4time_3dw,NX_L,NY_L
     1                          ,kmax,u_3d,v_3d,ext,istatus)
!                   call mean_wind(u_3d,v_3d,topo,NX_L,NY_L,NZ_L
!       1               ,dum1_array,dum2_array,dum3_array,idum1_array
!       1                                               ,umean,vmean,,,istatus)
                    grid_spacing_m = sqrt(
     1                 (  lat(1,2) - lat(1,1)                  )**2
     1               + ( (lon(1,2) - lon(1,1))*cosd(lat(1,1))  )**2
     1                                  )    * 111317. ! Grid spacing m

                    call advect(umean,vmean,radar_array
     1                  ,dummy_array,grid_spacing_m,NX_L,NY_L
     1                  ,radar_array_adv,float(laps_cycle_time),1.,lon)

                endif ! Pregenerated file


!               Display Advected Reflectivity Field
                call plot_cont(radar_array_adv,1e0,0.,chigh,cint_ref
     1         ,asc9_tim_r,'LAPS Max Reflectivity  1 HR Fcst',
     1          i_overlay,c_display,'nest7grid',lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

!           elseif(c_field .eq. 'nt')then
!               l_radar_read = .false.
!               goto2010

            endif ! c_field

        elseif( c_type .eq. 'sn')then
            mode = 1

            call make_fnam_lp(I4time_radar,asc9_tim_r,istatus)

!           call get__low__ref(grid_ra_ref,topo,NX_L,NY_L,NZ_L,radar_array)

!           Read in surface temp data
            var_2d = 'T'
            ext = 'lsx'
            call get_laps_2dgrid(i4time_radar,laps_cycle_time,i4time_tem
     1p,
     1      ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                          ,temp_2d,0,istatus)

            if(istatus .ne. 1)then
                write(6,*)' LAPS Sfc Temp not available'
                return
            endif

!           Read in surface dewpoint data
            var_2d = 'TD'
            ext = 'lsx'
            call get_laps_2d(i4time_temp,
     1  ext,var_2d,units_2d,comment_2d,
     1                          NX_L,NY_L,td_2d,istatus)

            if(istatus .ne. 1)then
                write(6,*)' LAPS Sfc Dewpoint not available'
                return
            endif

!           Read in surface pressure data
            var_2d = 'PS'
            ext = 'lsx'
            call get_laps_2d(i4time_temp,
     1  ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                          ,pres_2d,istatus)

            if(istatus .ne. 1)then
                write(6,*)' LAPS Sfc Pressure not available'
                return
            endif

            call get_tw_approx_2d(temp_2d,td_2d,pres_2d,NX_L,NY_L,tw_sfc
     1_k)

            call zs(radar_array,temp_2d,td_2d,pres_2d,tw_sfc_k,NX_L,NY_L
     1                                                  ,snow_2d)

!           c33_label = 'LAPS Snow Rate (liq equiv) in/hr '
!           scale = 1. / ((100./2.54) * 3600.) ! DENOM = (IN/HR) / (M/S)

            c33_label = 'LAPS Snowfall Rate         in/hr '
            scale = 1. / (10. * (100./2.54) * 3600.) ! DENOM = (IN snow/HR) / (M/S)

            call plot_cont(snow_2d,scale,
     1             0.,0.,-0.1,asc9_tim_r,c33_label,
     1          i_overlay,c_display,'nest7grid',lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

        elseif( c_type .eq. 'sa' .or. c_type .eq. 'pa' )then
            if(c_type .eq. 'sa')then
                write(6,1321)
1321            format('     ','Enter # of Hours of Snow Accumulation,',
     1          ' [-99 for Storm Total]     ','? '$)
                var_2d = 'STO'
            else
                write(6,1322)
1322            format('     ','Enter # of Hours of Precip Accumulation,
     1',
     1          ' [-99 for Storm Total]   ','? '$)
                var_2d = 'RTO'
            endif

            read(lun,*)r_hours

            write(6,*)

            ext = 'l1s'
            call get_directory(ext,directory,len_dir)

!           Cycle over at :28 after (if input time is not on the hour)
            if(i4time_ref .ne. (i4time_ref / 3600) * 3600)then
                i4time_ref1 = (i4time_ref-1680)/laps_cycle_time
     1                                        * laps_cycle_time
            else
                i4time_ref1 = i4time_ref
            endif

            if(r_hours .eq. -99.)then ! Storm Total
                write(6,*)' Getting Storm Total Accum from file'
                c9_string = 'Storm Tot'
                call get_laps_2dgrid(i4time_ref1,10000000,i4time_stm_tot
     1                  ,ext,var_2d
     1                  ,units_2d,comment_2d,NX_L,NY_L,accum_2d,0
     1                                                  ,istatus)
                write(6,*)' Storm Total was reset at ',comment_2d(1:9)
                call i4time_fname_lp(comment_2d(1:9),I4time_reset,istatu
     1s)
                istatus = 1
                num_hr_accum = (i4time_stm_tot - i4time_reset) / 3600
                i4time_accum = i4time_stm_tot
                i4time_end = i4time_stm_tot
                i4time_start = i4time_reset

!               encode(7,2017,c7_string)min(num_hr_accum,999)
                write(c7_string,2017)min(num_hr_accum,999)
2017            format(i4,' Hr')
                if(c_type .eq. 'sa')then
                    c33_label = 'LAPS Stm Tot Snow Acc (in)'//c7_string
                else
                    c33_label = 'LAPS Stm Tot Prcp Acc (in)'//c7_string
                endif

!           elseif(i4time_end .lt. 970677700)then
!               write(6,*)
!       1      ' Getting Entire Time Span of Accumulation from Radar Data, etc.'
!               encode(9,2029,c9_string)r_hours
!2029           format(f5.1,' Hr ')
!               call get_snow_accum(i4time_start,i4time_end,NX_L,NY_L,NZ_L
!       1               ,lat,lon,topo,grid_ra_ref,grid_ra_vel
!       1                            ,snow_2d,istatus)

            else ! Near Realtime - look for snow accumulation files
                if(i4time_now_gg() - i4time_ref1 .lt. 300)then ! Real Time Radar
                   !Find latest time of radar data
                    if('nest7grid' .ne. 'STORMFEST')then ! Read MHR packed data
                        c_filespec = c_filespec_ra
                    else
                        c_filespec = c_filespec_src
                    endif

                    call get_file_time(c_filespec,i4time_ref1,i4time_acc
     1um)
                else
!                   i4time_accum = (i4time_ref1+60) / 120 * 120 ! Rounded off time
                    i4time_accum = i4time_ref1
                endif

                i4time_end = i4time_accum
                i4time_interval = nint(r_hours * 3600.)
                i4time_start = i4time_end - i4time_interval

!               Round down to nearest cycle time
                i4time_endfile   = i4time_end  /laps_cycle_time*laps_cyc
     1le_time

!               Round up to nearest cycle time
                i4time_startfile = i4time_start/laps_cycle_time*laps_cyc
     1le_time

                if(i4time_start .gt. i4time_startfile)
     1          i4time_startfile = i4time_startfile + laps_cycle_time

                if(i4time_startfile .lt. i4time_endfile)then
                    istatus_file = 1
                    write(6,*)
     1     ' Looking for Storm Total Accumulations Stored In Files'

                    call get_laps_2d(i4time_endfile,ext,var_2d
     1          ,units_2d,comment_2d,NX_L,NY_L,accum_2d_buf,istatus_file
     1)
                    if(istatus_file .ne. 1)goto2100
                    comment_b = comment_2d(1:9)

                    call i4time_fname_lp(comment_2d(1:9),I4time_reset,is
     1tatus)
                    istatus = 1

                    if(i4time_startfile .ne. i4time_reset)then
                        call get_laps_2d(i4time_startfile,ext,var_2d
     1             ,units_2d,comment_2d,NX_L,NY_L,accum_2d,istatus_file)
                        if(istatus_file .ne. 1)goto2100
                        comment_a = comment_2d(1:9)

                        if(comment_a .ne. comment_b)then
                            write(6,*)' Storm Total was reset at '
     1                                            ,comment_b
!                           write(6,*)' Storm Total Clock was reset ',co
!    1mment_a
!    1                                                   ,' ',comment_b
                            write(6,*)
     1               ' Cannot subtract storm totals to get accumulation'
                            istatus_file = 0
                            goto2100
                        endif

                    else ! Reset time = Start File Time
                        write(6,*)
     1           ' Start File Time = Reset time, Set Init Accum to 0'
                        call zero(accum_2d,NX_L,NY_L)

                    endif

                    write(6,*)
     1     ' Subtracting Storm Totals to yield Accumulation over period'
                    call diff(accum_2d_buf,accum_2d,accum_2d,NX_L,NY_L)

                    do j = 1,NY_L
                    do i = 1,NX_L
                        if(accum_2d(i,j) .lt. 0.)then
                            write(6,*)' This should never happen:'
                            write(6,*)' Negative accum; Storm Total was
     1reset'
                            istatus_file = 0
                            goto2100
                        endif
                    enddo ! i
                    enddo ! j

                else
                    istatus_file = 0

                endif

2100            if(istatus_file .eq. 1)then ! Fill in ends of Pd with radar etc. data
                    if(i4time_start .lt. i4time_startfile)then
                        write(6,*)' Sorry, no L1S files present'
                        goto 1200

                    endif

                    if(i4time_end .gt. i4time_endfile)then
                        write(6,*)' Sorry, no L1S files present'
                        goto 1200

                    endif

                else ! Get entire time span from radar etc. data

                 write(6,*)' Sorry, no L1S files present'
                 goto 1200

                 write(6,*)
     1       ' Getting Entire Time Span of Accumulation from Radar Data,
     1 etc.'

                    if(c_type .eq. 'sa')then
                        call move(snow_2d,accum_2d,NX_L,NY_L)
                    else
                        call move(precip_2d,accum_2d,NX_L,NY_L)
                    endif

                endif

!               encode(9,2029,c9_string)r_hours
                write(c9_string,2029)r_hours
2029            format(f5.1,' Hr ')

                if(c_type .eq. 'sa')then
                    c33_label = 'LAPS '//c9_string//' Snow Accum  (in)
     1'
                else
                    c33_label = 'LAPS '//c9_string//' Prcp Accum  (in)
     1'
                endif

            endif

            if(istatus .ne. 1)goto1200

            call make_fnam_lp(I4time_accum,asc9_tim_r,istatus)

            scale = 1. / ((100./2.54)) ! DENOM = (IN/M)


            if(c_type .eq. 'pa')then
                if(abs(r_hours) .gt. 1.0)then
                    cint = -0.05
                else
                    cint = -0.01
                endif
                chigh = 50.
            else ! 'sa'
                if(abs(r_hours) .gt. 1.0)then
                    cint = -0.2
                else
                    cint = -0.1
                endif
                chigh = 200.
            endif

!           Eliminate "minor" maxima
            do j = 1,NY_L
            do i = 1,NX_L
                if(c_type .eq. 'pa')then
                    if(accum_2d(i,j)/scale .lt. 0.005)then
                        accum_2d(i,j) = 0.0
                    endif
!                   if(accum_2d(i,j)/scale .le. 0.0001)then
!                       accum_2d(i,j) = -1e-6
!                   endif
                else
                    if(accum_2d(i,j)/scale .lt. 0.05)then
                        accum_2d(i,j) = 0.0
                    endif
!                   if(accum_2d(i,j)/scale .le. 0.0001)then
!                       accum_2d(i,j) = -1e-6
!                   endif
                endif
            enddo ! i
            enddo ! j

            call plot_cont(accum_2d,scale,
     1             0.,chigh,cint,asc9_tim_r,c33_label,
     1          i_overlay,c_display,'nest7grid',lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

        elseif( c_type .eq. 'rx')then
            write(6,1311)
1311        format('     ','Enter # of Hours of Radar Data,',
     1          ' [-99 for Storm Total]     ','? '$)
            read(lun,*)r_hours

            write(6,*)

            i4time_accum = (i4time_ref+60) / 120 * 120 ! Rounded off time
            i4time_end = i4time_accum
            i4time_interval = nint(r_hours * 3600.)
            i4time_start = i4time_end - i4time_interval

            if(.true.)then ! Get entire time span from radar etc. data
                 write(6,*)
     1       ' Getting Entire Time Span of Accumulation from Radar Data,
     1 etc.'
                 call get_radar_max_pd(i4time_start,i4time_end
     1    ,NX_L,NY_L,NZ_L,lat,lon,topo,grid_ra_ref
     1        ,dummy_array,radar_array,frac_sum,istatus)

            endif

!           encode(9,2029,c9_string)r_hours
            write(c9_string,2029)r_hours
!2029        format(f5.1,' Hr ')

            c33_label = 'LAPS '//c9_string//' Reflctvty History '

            if(istatus .ne. 1)goto1200

            call make_fnam_lp(I4time_accum,asc9_tim_r,istatus)

            scale = 1.

            if(.false.)then
                write(6,*)' writing LRX field in current directory'
                directory = '[]'
                ext = 'lrx'
                var_2d = 'LRX'
                call put_laps_2d(i4time_accum,ext,var_2d
     1          ,units_2d,comment_2d,NX_L,NY_L,radar_array,istatus)
            endif

            call plot_cont(radar_array,scale,
     1             0.,80.,cint_ref,asc9_tim_r,c33_label,
     1          i_overlay,c_display,'nest7grid',lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

        elseif(c_type .eq. 't' .or. c_type .eq. 'pt')then
            write(6,1513)
1513        format('     Enter Level in mb',48x,'? '$)
            read(lun,*)k_level

!           if(istatus .ne. 1)goto1200

            if(c_type .eq. 'pt')then
                iflag_temp = 0 ! Returns Potential Temperature

                if(k_level .gt. 0)then
                    k_level = nint(zcoord_of_pressure(float(k_level*100)
     1))
                endif

                call get_temp_3d(i4time_ref,i4time_nearest,iflag_temp
     1                          ,NX_L,NY_L,NZ_L,temp_3d,istatus)

                call mklabel33(k_level,'  Potential Temp  K',c33_label)
                do i = 1,NX_L
                do j = 1,NY_L
                    temp_2d(i,j) = temp_3d(i,j,k_level)
                enddo ! j
                enddo ! i

            elseif(c_type .eq. 't ')then

                call get_temp_2d(i4time_ref,7200,i4time_nearest
     1                          ,k_level,NX_L,NY_L,temp_2d,istatus)

                do i = 1,NX_L
                do j = 1,NY_L
                    temp_2d(i,j) = temp_2d(i,j) - 273.15
                enddo ! j
                enddo ! i

                if(k_level .gt. 0)then
                    k_level = nint(zcoord_of_pressure(float(k_level*100)
     1))
                endif

                call mklabel33(k_level,'  Temperature     C',c33_label)

            endif

            call make_fnam_lp(i4time_nearest,asc9_tim_t,istatus)

            if(pressure_of_level(k_level) .le. 80000.)then
                clow =  0.
                chigh = 0.
                cint = 2.
            else
                clow =  0.
                chigh = 0.
                cint = 5.
            endif

            call plot_cont(temp_2d,1e0,clow,chigh,cint,
     1  asc9_tim_t,c33_label,
     1  i_overlay,c_display,'nest7grid',lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

        elseif(c_type .eq. 'hh')then
            write(6,1515)
1515        format('     Enter Temperature surface to display '
     1                      ,'height of (deg C)',11x,'? '$)
            read(lun,*)temp_in_c
            temp_in_k = temp_in_c + 273.15

            iflag_temp = 1 ! Returns Ambient Temperature

            call get_temp_3d(i4time_ref,i4time_nearest,iflag_temp
     1                          ,NX_L,NY_L,NZ_L,temp_3d,istatus)
!           if(istatus .ne. 1)goto1200

            write(c33_label,1516)nint(temp_in_c)
1516        format('LAPS Height of ',i3,'C Lvl (hft MSL)')

            do j = 1,NY_L
            do i = 1,NX_L
                height_2d(i,j) = 0.

                do k = 1,NZ_L-1
                    if(temp_3d(i,j,k  ) .gt. temp_in_k  .and.
     1               temp_3d(i,j,k+1) .le. temp_in_k       )then

!                       Desired Temperature occurs in this layer
                        frac_k = (       temp_in_k - temp_3d(i,j,k))/
     1                         (temp_3d(i,j,k+1) - temp_3d(i,j,k))
                        height_2d(i,j) =
     1             height_of_level(k)   * (1.0 - frac_k)
     1                 + height_of_level(k+1) * frac_k

                        height_2d(i,j) = height_2d(i,j) * 3.281

                    endif
                enddo ! k
            enddo ! i
            enddo ! j

            clow = 0.
            chigh = 0.
            cint = 5.

            call make_fnam_lp(i4time_nearest,asc9_tim_t,istatus)

            call plot_cont(height_2d,1e2,clow,chigh,cint,asc9_tim_t,
     1  c33_label,i_overlay,c_display,'nest7grid',lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

        elseif(c_type .eq. 'la' .or. c_type .eq. 'lj' .or.
     1       c_type .eq. 'sj' .or. c_type .eq. 'ls' .or.
     1       c_type .eq. 'ss' .or. c_type .eq. 'ci'                   )t
     1hen
            write(6,1514)
1514        format('     Enter Level in mb; OR [-1] for max in column'
     1                          ,21x,'? '$)
            read(lun,*)k_level
            k_mb = k_level

            if(k_level .gt. 0)then
               k_level = nint(zcoord_of_pressure(float(k_level*100)))
            endif

            i4time_lwc = i4time_ref/laps_cycle_time * laps_cycle_time

!           if(i4time_now_gg() - i4time_lwc .lt. 43200
!       1                               .and. c_type .eq. 'ls')then
                l_pregen = lapsplot_pregen
!           else
!               l_pregen = .false.
!           endif

            if(c_type .eq. 'la')then
                iflag_slwc = 1 ! Returns Adiabatic LWC
                if(k_level .gt. 0)then
                    call mklabel33(k_level,
     1          ' Adiabt LWC  g/m^3 ',c33_label)
                else
                    c33_label = 'LAPS Maximum Adiabatic LWC g/m^3 '
                endif

            elseif(c_type .eq. 'lj')then
                iflag_slwc = 2 ! Returns Adjusted LWC
                if(k_level .gt. 0)then
                    call mklabel33(k_level,
     1          ' Adjstd LWC  g/m^3 ',c33_label)
                else
                    c33_label = 'LAPS Maximum Adjusted  LWC g/m^3 '
                endif

            elseif(c_type .eq. 'sj')then
                iflag_slwc = 3 ! Returns Adjusted SLWC
                if(k_level .gt. 0)then
                    call mklabel33(k_level,
     1          ' Adjstd SLWC g/m^3 ',c33_label)
                else
                    c33_label = 'LAPS Maximum Adjusted SLWC g/m^3 '
                endif

            elseif(c_type .eq. 'ls')then
                iflag_slwc = 13 ! Returns New Smith - Feddes LWC
                if(k_level .gt. 0)then
                    call mklabel33(k_level,
     1          ' Smt-Fed LWC g/m^3 ',c33_label)
                else
                    c33_label = 'LAPS Max Smith-Feddes  LWC g/m^3 '
                endif

            elseif(c_type .eq. 'ci')then
                iflag_slwc = 13 ! Returns Cloud Ice
                if(k_level .gt. 0)then
                    call mklabel33(k_level,
     1          ' Smt-Fed ICE g/m^3 ',c33_label)
                else
                    c33_label = 'LAPS Max Smith-Feddes  ICE g/m^3 '
                endif

            elseif(c_type .eq. 'ss')then
                iflag_slwc = 14 ! Returns Smith - Feddes SLWC
                if(k_level .gt. 0)then
                    call mklabel33(k_level,
     1          'Smt-Fed SLWC g/m^3 ',c33_label)
                else
                    c33_label = 'LAPS Max Smith-Feddes SLWC g/m^3 '
                endif

            endif


            if(l_pregen)then
                write(6,*)' Getting pregenerated LWC file'
                if(c_type .ne. 'ci')then
                    var_2d = 'LWC'
                else
                    var_2d = 'ICE'
                endif
                ext = 'lwc'
                call get_directory(ext,directory,len_dir)

                if(k_mb .eq. -1)then ! Get 3D Grid
                  if(c_type .ne. 'ci')then
                    call get_laps_3dgrid(i4time_ref,86400,i4time_cloud,
     1          NX_L,NY_L,NZ_L,ext,var_2d
     1                  ,units_2d,comment_2d,slwc_3d,istatus)
                  else
                    call get_laps_3dgrid(i4time_ref,86400,i4time_cloud,
     1          NX_L,NY_L,NZ_L,ext,var_2d
     1                  ,units_2d,comment_2d,cice_3d,istatus)
                  endif

                else ! Get 2D horizontal slice from 3D Grid
                  if(c_type .ne. 'ci')then
                    call get_laps_2dgrid(i4time_ref,86400,i4time_cloud,
     1                  ext,var_2d
     1            ,units_2d,comment_2d,NX_L,NY_L,field_2d,k_mb,istatus)
                  else
                    call get_laps_2dgrid(i4time_ref,86400,i4time_cloud,
     1                  ext,var_2d
     1            ,units_2d,comment_2d,NX_L,NY_L,cice_2d,k_mb,istatus)
                  endif

                endif

            else ! Calculate on the Fly

            endif ! L_pregen

            call make_fnam_lp(i4time_lwc,asc9_tim_t,istatus)

            clow = 0.
            chigh = 0.
            cint = -0.1

            if(k_level .gt. 0)then ! Plot SLWC on const pressure sfc
               if(c_type .ne. 'ci')then
                 if(l_pregen)then
                   call plot_cont(field_2d,1e-3,
     1               clow,chigh,cint,asc9_tim_t,c33_label,
     1          i_overlay,c_display,'nest7grid',lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)
                 else
                   call plot_cont(slwc_3d(1,1,k_level),1e-3,
     1               clow,chigh,cint,asc9_tim_t,c33_label,
     1          i_overlay,c_display,'nest7grid',lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)
                 endif
               else ! c_type .ne. 'ci'
                 if(l_pregen)then
                   call plot_cont(cice_2d,1e-3,
     1               clow,chigh,cint,asc9_tim_t,c33_label,
     1          i_overlay,c_display,'nest7grid',lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)
                 else
                   call plot_cont(cice_3d(1,1,k_level),1e-3,
     1               clow,chigh,cint,asc9_tim_t,c33_label,
     1          i_overlay,c_display,'nest7grid',lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)
                 endif
               endif ! c_type

            else ! Find Maximum value in column
               do j = 1,NY_L
               do i = 1,NX_L
                   column_max(i,j) = -1e-30
                   if(c_type .ne. 'ci')then
                     do k = 1,NZ_L
                       column_max(i,j) = max(column_max(i,j),slwc_3d(i,j
     1,k))
                     enddo ! k
                   else
                     do k = 1,NZ_L
                       column_max(i,j) = max(column_max(i,j),cice_3d(i,j
     1,k))
                     enddo ! k
                   endif
               enddo ! i
               enddo ! j
               call plot_cont(column_max,1e-3,
     1               clow,chigh,cint,asc9_tim_t,c33_label,
     1          i_overlay,c_display,'nest7grid',lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

            endif

        elseif(c_type .eq. 'mv' .or. c_type .eq. 'ic')then
            write(6,1514)

            read(lun,*)k_level
            k_mb = k_level

            if(k_level .gt. 0)then
               k_level = nint(zcoord_of_pressure(float(k_level*100)))
            endif

            i4time_lwc = i4time_ref/laps_cycle_time * laps_cycle_time

            if(c_type .eq. 'mv')then
                if(k_level .gt. 0)then
                    call mklabel33(k_level,'     MVD     m^-6  ',c33_lab
     1el)
                else
                    c33_label = 'LAPS Mean Volume Diameter  m^-6  '
                endif

                write(6,*)' Getting pregenerated LMD file'
                var_2d = 'LMD'
                ext = 'lmd'

            elseif(c_type .eq. 'ic')then
                if(k_level .gt. 0)then
                    call mklabel33(k_level,'   Icing Index     ',c33_lab
     1el)
                else
                    c33_label = '        LAPS Icing Index         '
                endif

                write(6,*)' Getting pregenerated LRP file'
                var_2d = 'LRP'
                ext = 'lrp'

            endif ! c_type .eq. 'ic'

            if(k_mb .eq. -1)then ! Get 3D Grid
                call get_laps_3dgrid(i4time_ref,10000000
     1                                  ,i4time_cloud
     1                                  ,NX_L,NY_L,NZ_L,ext,var_2d
     1                                  ,units_2d,comment_2d,slwc_3d
     1                                  ,istatus)

            else ! Get 2D horizontal slice from 3D Grid
                call get_laps_2dgrid(i4time_ref,10000000
     1                                  ,i4time_cloud
     1                                  ,ext,var_2d
     1                                  ,units_2d,comment_2d,NX_L
     1                                  ,NY_L,field_2d,k_mb,istatus)

            endif


            call make_fnam_lp(i4time_cloud,asc9_tim_t,istatus)

            if(c_type .eq. 'mv')then
                clow = 10.
                chigh = 26.
                cint = 2.

                if(k_level .gt. 0)then ! Plot MVD on const pressure sfc
                   if(.true.)then
                       call plot_cont(mvd_2d,0.9999e-6,
     1                   clow,chigh,cint,asc9_tim_t,c33_label,
     1                   i_overlay,c_display,'nest7grid',lat,lon,jdot,
     1                   NX_L,NY_L,r_missing_data,laps_cycle_time)
                   else
                       call plot_cont(mvd_3d(1,1,k_level),0.9999e-6,
     1                   clow,chigh,cint,asc9_tim_t,c33_label,
     1                   i_overlay,c_display,'nest7grid',lat,lon,jdot,
     1                   NX_L,NY_L,r_missing_data,laps_cycle_time)
                   endif

                else ! Find Maximum value in column
                   do j = 1,NY_L
                   do i = 1,NX_L
                       column_max(i,j) = -1e-30
                       do k = 1,NZ_L
                           column_max(i,j) = max(column_max(i,j)
     1                                          ,mvd_3d(i,j,k))
                       enddo ! k
                   enddo ! i
                   enddo ! j

                   call plot_cont(column_max,0.9999e-6,
     1               clow,chigh,cint,asc9_tim_t,c33_label,
     1               i_overlay,c_display,'nest7grid',lat,lon,jdot,
     1               NX_L,NY_L,r_missing_data,laps_cycle_time)

                endif

            elseif(c_type .eq. 'ic')then
                clow = 0.
                chigh = 10.
                cint = 1.0

                if(k_level .gt. 0)then ! Plot on const pressure sfc
                   if(.true.)then
                       call plot_cont(field_2d,1e0,clow,chigh,cint
     1                      ,asc9_tim_t
     1                      ,c33_label,i_overlay,c_display
     1                      ,'nest7grid',lat,lon,jdot
     1                      ,NX_L,NY_L,r_missing_data,laps_cycle_time)

                   endif

                else ! Find Maximum value in column
                   if(.true.)then
                       do j = 1,NY_L
                       do i = 1,NX_L
                           column_max(i,j) = -1e-30
                           do k = 1,NZ_L
                               if(slwc_3d(i,j,k) .gt. 0.)column_max(i,j)
     1                     = max(column_max(i,j),slwc_3d(i,j,k)+.01)
                           enddo ! k
                       enddo ! i
                       enddo ! j
                   endif

                   call plot_cont(column_max,1e0,
     1                  clow,chigh,cint,asc9_tim_t,c33_label,
     1                  i_overlay,c_display,'nest7grid',lat,lon,jdot,
     1                  NX_L,NY_L,r_missing_data,laps_cycle_time)

                endif ! k_level

            endif ! c_type

        elseif(c_type .eq. 'tc')then
1524        write(6,1517)
1517        format('     Enter Lvl (mb); OR [0] 2D cldtyp'
!    1          ,' [-1] low cloud,'
!    1          ,' [-2] high cloud'
     1          ,' ? '$)

1525        read(lun,*)k_level
            k_mb = k_level

            if(k_level .lt. 0)then
                write(6,*)' Try Again'
                goto1524
            endif

            if(.true.)then ! Read 2D cloud type field

                if(k_level .gt. 0)then ! Read from 3-D cloud type
                   k_level =
     1                   nint(zcoord_of_pressure(float(k_level*100)))
                    ext = 'lty'
                    var_2d = 'CTY'
                    call mklabel33
     1                    (k_level,'     Cloud Type    ',c33_label)

                else                   ! Read from 2-D cloud type
                    ext = 'lct'
                    var_2d = 'SCT'
                    c33_label = '      LAPS    2-D Cloud Type     '

                endif


                call get_laps_2dgrid(i4time_ref,laps_cycle_time*100
     1            ,i4time_cloud,ext,var_2d
     1            ,units_2d,comment_2d,NX_L,NY_L,field_2d,k_mb,istatus)

                IF(istatus .ne. 1 .and. istatus .ne. -1)THEN
                    write(6,*)' Error reading cloud type'
                    goto 1200
                endif

                call make_fnam_lp(i4time_cloud,asc9_tim,istatus)

!               Convert from real to byte
                do i = 1,NX_L
                do j = 1,NY_L
!                   Convert to byte
                    b_array(i,j) = i4_to_byte(int(field_2d(i,j)))
                enddo ! i
                enddo ! j


                call plot_cldpcp_type(b_array
     1                ,asc9_tim,c33_label,c_type,k,i_overlay,c_display
     1                ,lat,lon,idum1_array,'nest7grid'
     1                ,NX_L,NY_L,laps_cycle_time,jdot)

            else ! OLD ARCHAIC CODE

            endif ! k_level .eq. 0


        elseif(c_type .eq. 'tp')then
1624        write(6,1617)
1617        format('     Enter Level in mb; [0] for surface,'
     1          ,' OR [-1] for sfc thresholded: ','? '$)

1625        read(lun,*)k_level
            k_mb = k_level

            if(k_level .gt. 0)then
                call mklabel33(k_level,'    Precip Type    ',c33_label)
            elseif(k_level .eq.  0)then
                c33_label = 'LAPS Sfc Precip Type   (nothresh)'
            elseif(k_level .eq. -1)then
                c33_label = 'LAPS Sfc Precip Type   (thresh)  '
            endif

            if(k_level .eq. -1)then
                var_2d = 'PTT'
                k_level = 0
            else
                var_2d = 'PTY'
            endif

            if(k_level .gt. 0)then
               k_level = nint(zcoord_of_pressure(float(k_level*100)))
            endif

            i4time_pcp = i4time_ref/laps_cycle_time * laps_cycle_time

            l_precip_pregen = .true.

            if(k_level .gt. 0)then ! Plot Precip Type on const pressure sfc
                if(l_precip_pregen)then ! Read pregenerated field

                    write(6,*)' Reading pregenerated precip type field'
                    ext = 'lty'
                    call get_laps_2dgrid(i4time_pcp,laps_cycle_time,i4ti
     1me_nearest,
     1                  ext,var_2d
     1        ,units_2d,comment_2d,NX_L,NY_L,field_2d,k_mb,istatus)

!                   Convert from real to byte
                    do i = 1,NX_L
                    do j = 1,NY_L
!                       Convert to integer
                        iarg = int(field_2d(i,j)) * 16
!                       Convert to byte
                        pcp_type_2d(i,j) = i4_to_byte(iarg)
                    enddo ! i
                    enddo ! j

                    call make_fnam_lp(i4time_nearest,asc9_tim,istatus)

                    call plot_cldpcp_type(pcp_type_2d
     1         ,asc9_tim,c33_label,c_type,k_level,i_overlay,c_display
     1         ,lat,lon,idum1_array,'nest7grid'
     1     ,NX_L,NY_L,laps_cycle_time,jdot)

                else
                    call plot_cldpcp_type(cldpcp_type_3d(1,1,k_level)
     1         ,asc9_tim,c33_label,c_type,k_level,i_overlay,c_display
     1         ,lat,lon,idum1_array,'nest7grid'
     1     ,NX_L,NY_L,laps_cycle_time,jdot)

                endif

            elseif(k_level .eq. 0)then ! Extract Surface Precip Type Field

                if(l_precip_pregen)then

                  ! Read SFC precip type from lty field
                    write(6,*)' Reading pregenerated SFC precip type fie
     1ld ',var_2d
!                   var_2d was defined earlier in the if block
                    ext = 'lct'
                    call get_laps_2dgrid(i4time_pcp,laps_cycle_time,i4ti
     1me_temp,
     1          ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                          ,field_2d,0,istatus)

                    if(istatus .ne. 1)goto1200

                    call make_fnam_lp(i4time_temp,asc9_tim,istatus)

!                   Convert from real to byte
                    do i = 1,NX_L
                    do j = 1,NY_L
                        iarg = field_2d(i,j) * 16 ! Code into left 4 bits
                        pcp_type_2d(i,j) = i4_to_byte(iarg)
                    enddo ! i
                    enddo ! j

                else ! Calculate Precip Type on the Fly
!                   Read in surface temp data
                    var_2d = 'T'
                    ext = 'lsx'
                    call get_laps_2dgrid(i4time_pcp,laps_cycle_time,i4ti
     1me_temp,
     1          ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                          ,temp_2d,0,istatus)

                    if(istatus .ne. 1)then
                        write(6,*)' LAPS Sfc Temp not available'
                        goto1200
                    endif

!                   Read in surface dewpoint data
                    var_2d = 'TD'
                    ext = 'lsx'
                    call get_laps_2d(i4time_pcp,
     1              ext,var_2d,units_2d,comment_2d,
     1                          NX_L,NY_L,td_2d,istatus)

                    if(istatus .ne. 1)then
                        write(6,*)' LAPS Sfc Dewpoint not available'
                        goto1200
                    endif

                    call get_sfc_preciptype(pres_2d,temp_2d,td_2d,cldpcp
     1_type_3d
     1                                  ,pcp_type_2d,NX_L,NY_L,NZ_L)

!                   call make_fnam_lp(i4time_pcp,asc9_tim,istatus)

                endif ! l_precip_pregen

                call plot_cldpcp_type(pcp_type_2d
     1     ,asc9_tim,c33_label,c_type,k,i_overlay,c_display
     1     ,lat,lon,idum1_array,'nest7grid'
     1     ,NX_L,NY_L,laps_cycle_time,jdot)

            endif ! k_level


!       elseif(c_type .eq. 'cl')then
!           iflag_temp = 1 ! Returns Ambient Temp (K)
!            call get_temp_3d(i4time_ref,i4time_nearest,iflag_temp
!       1                               ,NX_L,NY_L,NZ_L,temp_3d,istatus)
!           if(istatus .ne. 1)goto1200

!           call make_fnam_lp(i4time_nearest,asc9_tim_t,istatus)

!           write(6,*)' Site?'

!           read(lun,701)c3_site
!701         format(a3)

!           call compare_classdat(temp_3d,asc9_tim_t,i4time_ref,c3_site)

        elseif(c_type .eq. 'ia' .or. c_type .eq. 'ij'
     1                        .or. c_type .eq. 'is')then

          if(.false.)then

          else
            ext = 'lil'
            var_2d = 'LIL'
            call get_laps_2dgrid(i4time_ref,86400,i4time_cloud,
     1                  ext,var_2d
     1            ,units_2d,comment_2d,NX_L,NY_L,column_max,0,istatus)

            call make_fnam_lp(i4time_cloud,asc9_tim_t,istatus)
            c33_label = 'LAPS Integrated Smith-Fed LWC mm '

          endif ! false

          clow = 0.
          chigh = +0.
          cint = -0.1

          call plot_cont(column_max,1e3,
     1          clow,chigh,cint,asc9_tim_t,c33_label,
     1                  i_overlay,c_display,'nest7grid',lat,lon,jdo
     1t,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

        elseif(c_type .eq. 'pe' .or. c_type .eq. 'ne')then
          if(lapsplot_pregen)then ! Pregenerated version
            ext = 'lsx'

            if(c_type .eq. 'pe')then
                var_2d = 'PBE'

                call get_laps_2dgrid(i4time_ref,10000000,i4time_temp,
     1          ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                          ,pbe_2d,0,istatus)

            else
                var_2d = 'NBE'

                call get_laps_2dgrid(i4time_ref,10000000,i4time_temp,
     1          ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                          ,nbe_2d,0,istatus)

            endif

            call make_fnam_lp(i4time_temp,asc9_tim_t,istatus)

          else
            iflag_temp = 1 ! Returns Ambient Temperature

            call get_temp_3d(i4time_ref,i4time_temp,iflag_temp
     1                          ,NX_L,NY_L,NZ_L,temp_3d,istatus)
!           if(istatus .ne. 1)goto1200

!           Read in surface temp data
            var_2d = 'T'
            ext = 'lsx'
            call get_laps_2d(i4time_temp,
     1  ext,var_2d,units_2d,comment_2d,
     1                          NX_L,NY_L,temp_2d,istatus)

            if(istatus .ne. 1)then
                write(6,*)' LAPS Sfc Temp not available'
                return
            endif

!           Read in surface dewpoint data
            var_2d = 'TD'
            ext = 'lsx'
            call get_laps_2d(i4time_temp,
     1  ext,var_2d,units_2d,comment_2d,
     1                          NX_L,NY_L,td_2d,istatus)

            if(istatus .ne. 1)then
                write(6,*)' LAPS Sfc Dewpoint not available'
                goto1200
            endif

!           Read in SFC pressure data
            var_2d = 'PS'
            ext = 'lsx'
            call get_laps_2d(i4time_temp,
     1  ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                          ,pres_2d,istatus)
            if(istatus .ne. 1)then
                write(6,*)' LAPS SFC pressure not available'
                goto1200
            endif

            do k = 1,NZ_L ! Calculate pressure at each level
                p_1d_pa(k) = pressure_of_level(k) ! Pressure at each level
            enddo ! k

!           Calculate a 3-D Height Field
            call get_heights_hydrostatic(temp_3d,pres_2d,topo
     1          ,dum1_array,dum2_array,dum3_array,dum4_array
     1          ,NX_L,NY_L,NZ_L
     1                                          ,heights_3d)

!           Get PBE and NBE - Make sure t_sfc_k(i,j) >= td_sfc_k(i,j)
            call laps_be(NX_L,NY_L,NZ_L
     1  ,temp_2d,td_2d,pres_2d,temp_3d,heights_3d,p_1d_pa,topo,pbe_2d,nb
     1e_2d)

          endif

          scale = 1.
          call make_fnam_lp(i4time_temp,asc9_tim_t,istatus)

          if(c_type .eq. 'pe')then
              c33_label = 'LAPS CAPE                (J/KG)  '
              clow = 0.
              chigh = 8000.
              cint = +400.
              call plot_cont(pbe_2d,scale,clow,chigh,cint,asc9_tim_t,c33
     1_label,
     1                  i_overlay,c_display,'nest7grid',lat,lon,jdo
     1t,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

          elseif(c_type .eq. 'ne')then
!             Change flag value (for now)
              do i = 1,NX_L
              do j = 1,NY_L
                  if(nbe_2d(i,j) .eq. -1e6)nbe_2d(i,j) = r_missing_data
                  if(nbe_2d(i,j) .ge. -2. .and. nbe_2d(i,j) .ne. r_missi
     1ng_data)
     1                                   nbe_2d(i,j) = +.0001
              enddo ! j
              enddo ! i

              c33_label = 'LAPS CIN                 (J/KG)  '
              clow = -500 !   0.
              chigh = 0.  !   0.
              cint = 50.  ! -10.
              call plot_cont(nbe_2d,scale,clow,chigh,cint,asc9_tim_t,c33
     1_label,
     1                  i_overlay,c_display,'nest7grid',lat,lon,jdo
     1t,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)
          endif

        elseif(c_type .eq. 'sh')then
            write(6,1513)
            read(lun,*)k_level
            if(k_level .gt. 0)then
               k_level = nint(zcoord_of_pressure(float(k_level*100)))
            endif

            var_2d = 'SH '
            ext = 'lq3'

            call get_laps_3dgrid
     1  (i4time_ref,1000000,i4time_nearest,NX_L,NY_L,NZ_L
     1          ,ext,var_2d,units_2d,comment_2d
     1                                  ,q_3d,istatus)

            call mklabel33(k_level,' Spec Humidity x1e3',c33_label)

            clow = 0.
            chigh = +40.
            cint = 0.2
            cint = -1.

            call make_fnam_lp(i4time_nearest,asc9_tim_t,istatus)

!           do i = 1,NX_L
!           do j = 1,NY_L      
!               if(q_3d(i,j,k_level) .eq. r_missing_data)then
!                   q_3d(i,j,k_level) = 0.
!               endif              
!           enddo ! j
!           enddo ! i

            call plot_cont(q_3d(1,1,k_level),1e-3,clow,chigh,cint,
     1  asc9_tim_t,c33_label,i_overlay,c_display,'nest7grid'
     1          ,lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

        elseif(c_type .eq. 'rh')then
            write(6,1513)
            read(lun,*)k_level
            if(k_level .gt. 0)then
               k_level = nint(zcoord_of_pressure(float(k_level*100)))
            endif

            var_2d = 'RHL'
            ext = 'lh3'

            call get_laps_3dgrid
     1  (i4time_ref,1000000,i4time_nearest,NX_L,NY_L,NZ_L
     1          ,ext,var_2d,units_2d,comment_2d
     1                                  ,rh_3d,istatus)

            call mklabel33(k_level,' Relative Humidity %',c33_label)

            clow = 0.
            chigh = +100.
            cint = 10.

            call make_fnam_lp(i4time_nearest,asc9_tim_t,istatus)

!           do i = 1,NX_L
!           do j = 1,NY_L      
!               if(rh_3d(i,j,k_level) .eq. r_missing_data)then
!                  rh_3d(i,j,k_level) = 0.
!               endif              
!           enddo ! j
!           enddo ! i

            call plot_cont(rh_3d(1,1,k_level),1e0,clow,chigh,cint,
     1  asc9_tim_t,c33_label,i_overlay
     1          ,c_display,'nest7grid',lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

        elseif(c_type .eq. 'hy')then
            write(6,1513)
            read(lun,*)k_level
            if(k_level .gt. 0)then
               k_level = nint(zcoord_of_pressure(float(k_level*100)))
            endif

            iflag_temp = 1 ! Returns Ambient Temperature
            call get_temp_3d(i4time_ref,i4time_nearest,iflag_temp
     1                          ,NX_L,NY_L,NZ_L,temp_3d,istatus)

!           Read in SFC pressure
            i4time_tol = 0
            var_2d = 'PS'
            ext = 'lsx'
            call get_laps_2dgrid(i4time_nearest,i4time_tol,i4time_neares
     1t,
     1          ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                  ,pres_2d,0,istatus)
            IF(istatus .ne. 1)THEN
                write(6,*)' Error Reading Surface Pres Analyses'
     1                  ,' - no hydrostatic heights calculated'
                goto1200
            endif

            call get_heights_hydrostatic(temp_3d,pres_2d,topo,
     1          dum1_array,dum2_array,dum3_array,dum4_array,
     1                                  NX_L,NY_L,NZ_L,heights_3d)

!           call get_laps_heights(i4time_ref,temp_3d,heights_3d)

            call mklabel33(k_level,' LAPS Heights    dm',c33_label)

            clow = 0.
            chigh = 0.
            cint = 1. ! 3.

            i4time_heights = i4time_nearest

            call make_fnam_lp(i4time_heights,asc9_tim_t,istatus)

            call plot_cont(heights_3d(1,1,k_level),1e1,clow,chigh,cint,
     1  asc9_tim_t,c33_label,i_overlay,c_display,'nest7grid'
     1          ,lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

        elseif(c_type .eq. 'hb')then

            var_2d = 'HT'

            call make_fnam_lp(i4time_ref,asc9_tim_t,istatus)

            ext = 'lga'
            call get_directory(ext,directory,len_dir)

            write(6,*)' Enter yydddhhmmHHMM for lga file'
            read(5,211)a13_time
 211        format(a13)
            call get_fcst_times(a13_time,I4TIME,i4_valid,i4_fn)

            write(6,1513)
            read(lun,*)k_mb
            k_level = nint(zcoord_of_pressure(float(k_mb*100)))
            k_mb    = nint(pressure_of_level(k_level) / 100.)

            CALL READ_LAPS(I4TIME,i4_valid,DIRECTORY,EXT,NX_L,NY_L,1,1,       
     1          VAR_2d,k_mb,LVL_COORD_2d,UNITS_2d,COMMENT_2d,
     1          field_2d,ISTATUS)

!           call get_laps_2dgrid(i4time_ref,laps_cycle_time*100,
!    1              i4time_heights,
!    1              ext,var_2d,units_2d,comment_2d,NX_L,NY_L
!    1                                     ,field_2d,k_mb,istatus)

            IF(istatus .ne. 1)THEN
                write(6,*)' Error Reading Background Height Analysis'
                goto1200
            endif

            call mklabel33(k_level,' LAPS Heights    dm',c33_label)

            clow = 0.
            chigh = 0.
            cint = 1. ! 3.

            call make_fnam_lp(i4time,asc9_tim_t,istatus)

            call plot_cont(field_2d,1e1,clow,chigh,cint,
     1  asc9_tim_t,c33_label,i_overlay,c_display,'nest7grid'
     1                  ,lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

        elseif(c_type .eq. 'ht')then
            write(6,1513)
            read(lun,*)k_mb

            k_level = nint(zcoord_of_pressure(float(k_mb*100)))

            var_2d = 'HT'

            call make_fnam_lp(i4time_ref,asc9_tim_t,istatus)

            ext = 'lt1'

            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100,i4time_h
     1eights,
     1              ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                     ,field_2d,k_mb,istatus)

            IF(istatus .ne. 1)THEN
                write(6,*)' Error Reading LAPS Height Analysis'
                goto1200
            endif

            call mklabel33(k_level,' LAPS Heights    dm',c33_label)

            clow = 0.
            chigh = 0.
            cint = 1. ! 3.

            call make_fnam_lp(i4time_heights,asc9_tim_t,istatus)

            call plot_cont(field_2d,1e1,clow,chigh,cint,
     1  asc9_tim_t,c33_label,i_overlay,c_display,'nest7grid'
     1                                          ,lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

        elseif(c_type .eq. 'pw')then
            var_2d = 'TPW'
            ext = 'lh4'
            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100,i4time_p
     1w,
     1              ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                     ,field_2d,0,istatus)

            IF(istatus .ne. 1)THEN
                write(6,*)' Error Reading Precipitable Water'
                goto1200
            endif

            c33_label = 'LAPS Total Precipitable Water  cm'

            clow = 0.
            chigh = 15.
            cint = .25

            call make_fnam_lp(i4time_pw,asc9_tim_t,istatus)

            call plot_cont(field_2d,1e-2,clow,chigh,cint,
     1  asc9_tim_t,c33_label,i_overlay,c_display,'nest7grid'
     1                                          ,lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

        elseif(c_type .eq. 'tt')then
            var_2d = 'T'
            ext = 'lsx'
            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100,i4time_p
     1w,
     1              ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                     ,field_2d,0,istatus)

!           K to F
            do i = 1,NX_L
            do j = 1,NY_L
                field_2d(i,j) = ((field_2d(i,j)-273.15) * 1.8) + 32.
            enddo ! j
            enddo ! i

            IF(istatus .ne. 1)THEN
                write(6,*)' Error Reading Surface Temps'
                goto1200
            endif

            c33_label = 'LAPS Sfc Temperature     (F)     '

            clow = -50.
            chigh = +120.
            cint = 5.

            call make_fnam_lp(i4time_pw,asc9_tim_t,istatus)

            call plot_cont(field_2d,1e-0,clow,chigh,cint,
     1  asc9_tim_t,c33_label,i_overlay,c_display,'nest7grid'
     1                                          ,lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

        elseif(c_type .eq. 'hi')then
            var_2d = 'HI'
            ext = 'lsx'
            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100,i4time_p
     1w,
     1              ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                     ,field_2d,0,istatus)

!           K to F
            do i = 1,NX_L
            do j = 1,NY_L
                field_2d(i,j) = ((field_2d(i,j)-273.15) * 1.8) + 32.
            enddo ! j
            enddo ! i

            IF(istatus .ne. 1)THEN
                write(6,*)' Error Reading Heat Index'
                goto1200
            endif

            c33_label = 'LAPS Heat Index          (F)     '

            clow = +50.
            chigh = +150.
            cint = 5.

            call make_fnam_lp(i4time_pw,asc9_tim_t,istatus)

            call plot_cont(field_2d,1e-0,clow,chigh,cint,
     1  asc9_tim_t,c33_label,i_overlay,c_display,'nest7grid'
     1                                          ,lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

        elseif(c_type .eq. 'td')then
            var_2d = 'TD'
            ext = 'lsx'
            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100,i4time_p
     1w,
     1              ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                     ,field_2d,0,istatus)

!           K to F
            do i = 1,NX_L
            do j = 1,NY_L
                field_2d(i,j) = ((field_2d(i,j)-273.15) * 1.8) + 32.
            enddo ! j
            enddo ! i

            IF(istatus .ne. 1)THEN
                write(6,*)' Error Reading Surface Td'
                goto1200
            endif

            c33_label = 'LAPS Sfc Dew Point       (F)     '

            clow = -50.
            chigh = +120.
            cint = 5.

            call make_fnam_lp(i4time_pw,asc9_tim_t,istatus)

            call plot_cont(field_2d,1e-0,clow,chigh,cint,
     1  asc9_tim_t,c33_label,i_overlay,c_display,'nest7grid'
     1                                          ,lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

        elseif(c_type .eq. 'mc')then
            var_2d = 'MRC'
            ext = 'lsx'
            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100,i4time_p
     1w,
     1              ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                     ,field_2d,0,istatus)

            IF(istatus .ne. 1)THEN
                write(6,*)' Error Reading Surface Moisture Convergence'
                goto1200
            endif

            c33_label = 'LAPS Sfc Mstr Flux Conv  (x 1e-4)'

            clow = -100.
            chigh = +100.
            cint = 5.

            call make_fnam_lp(i4time_pw,asc9_tim_t,istatus)

            call plot_cont(field_2d,1e-4,clow,chigh,cint,
     1  asc9_tim_t,c33_label,i_overlay,c_display,'nest7grid'
     1                                          ,lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)
        elseif(c_type .eq. 'ws')then ! surface wind
            ext = 'lsx'
            var_2d = 'U'
            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100,i4time_p
     1w,
     1              ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                     ,u_2d,0,istatus)
            var_2d = 'V'
            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100,i4time_p
     1w,
     1              ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                     ,v_2d,0,istatus)

            IF(istatus .ne. 1)THEN
                write(6,*)' Error Reading Surface Wind'
                goto1200
            endif

            c33_label = 'LAPS Surface Wind            (kt)'

            if(max(NX_L,NY_L) .gt. 50)then
                interval = 2
            else
                interval = 1
            endif

            size = float(interval) * .15

            call make_fnam_lp(i4time_pw,asc9_tim_t,istatus)

!           Rotate sfc winds from grid north to true north
            do i = 1,NX_L
            do j = 1,NY_L
                u_grid = u_2d(i,j)
                v_grid = v_2d(i,j)
                call uvgrid_to_uvtrue(  u_grid,
     1                                  v_grid,
     1                                  u_true,
     1                                  v_true,
     1                                  lon(i,j) )
                u_2d(i,j) = u_true
                v_2d(i,j) = v_true
            enddo ! j
            enddo ! i

            call plot_barbs(u_2d,v_2d,lat,lon,topo,size,interval
     1                     ,asc9_tim_t,c33_label,c_field,k_level
     1                     ,i_overlay,c_display,'nest7grid'
     1                     ,NX_L,NY_L,NZ_L,grid_ra_ref,grid_ra_vel
     1                     ,NX_L,NY_L,r_missing_data,laps_cycle_time
     1                     ,jdot)

        elseif(c_type .eq. 'p' .or. c_type .eq. 'pm')then ! 1500m or MSL Pres
            if(c_type .eq. 'p')then
                var_2d = 'P'
                c33_label = 'LAPS 1500m Pressure          (mb)'
            else
                var_2d = 'MSL'
                c33_label = 'LAPS MSL Pressure            (mb)'
            endif

            ext = 'lsx'
            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100
     1                          ,i4time_pw,ext,var_2d,units_2d
     1                          ,comment_2d,NX_L,NY_L,field_2d,0
     1                          ,istatus)

            IF(istatus .ne. 1)THEN
                write(6,*)' Error Reading Surface ',var_2d
                goto1200
            endif

            pres_low_mb = 2000.
            pres_high_mb = 0.

!           Get range of pressures (in mb)
            do i = 1,NX_L
            do j = 1,NY_L
                pres_low_mb = min(pres_low_mb,field_2d(i,j)/100.)
                pres_high_mb = max(pres_high_mb,field_2d(i,j)/100.)
            enddo ! j
            enddo ! i

            cint = 2.
            icint = cint
            clow = (int(pres_low_mb) / icint) * icint
            chigh = (int(pres_high_mb) / icint) * icint + icint

            call make_fnam_lp(i4time_pw,asc9_tim_t,istatus)

            call plot_cont(field_2d,1e+2,clow,chigh,cint,asc9_tim_t
     1                    ,c33_label,i_overlay,c_display,'nest7grid'
     1                    ,lat,lon,jdot
     1                    ,NX_L,NY_L,r_missing_data,laps_cycle_time)

        elseif(c_type .eq. 'ps')then ! Surface Pressure
            var_2d = 'PS'
            ext = 'lsx'
            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100,i4time_p
     1w,
     1              ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                     ,field_2d,0,istatus)

            IF(istatus .ne. 1)THEN
                write(6,*)' Error Reading Surface ',var_2d
                goto1200
            endif

            c33_label = 'LAPS Surface Pressure        (mb)'

            clow = 650.
            chigh = +1100.
            cint = 5.

            call make_fnam_lp(i4time_pw,asc9_tim_t,istatus)

            call plot_cont(field_2d,1e+2,clow,chigh,cint,
     1  asc9_tim_t,c33_label,i_overlay,c_display,'nest7grid'
     1                                          ,lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

        elseif(c_type .eq. 'vv')then
            var_2d = 'VV'
            ext = 'lsx'
            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100,i4time_p
     1w,
     1              ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                     ,field_2d,0,istatus)

            IF(istatus .ne. 1)THEN
                write(6,*)' Error Reading Surface ',var_2d
                goto1200
            endif

            c33_label = 'LAPS Sfc Vert Velocity     (cm/s)'

            clow = -200.
            chigh = +200.
            cint = 10.

            call make_fnam_lp(i4time_pw,asc9_tim_t,istatus)

            call plot_cont(field_2d,1e-2,clow,chigh,cint,
     1  asc9_tim_t,c33_label,i_overlay,c_display,'nest7grid'
     1                                          ,lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)
        elseif(c_type .eq. 'hu')then
            var_2d = 'RH'
            ext = 'lsx'
            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100,i4time_p
     1w,
     1              ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                     ,field_2d,0,istatus)

            IF(istatus .ne. 1)THEN
                write(6,*)' Error Reading Surface ',var_2d
                goto1200
            endif

!           c33_label = 'LAPS Sfc  Rel Hum       (PERCENT)'
            c33_label = 'LAPS Sfc Rel Humidity   (PERCENT)'

            clow = 0.
            chigh = +100.
            cint = 10.

            call make_fnam_lp(i4time_pw,asc9_tim_t,istatus)

            call plot_cont(field_2d,1e-0,clow,chigh,cint,
     1  asc9_tim_t,c33_label,i_overlay,c_display,'nest7grid'
     1                                          ,lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

        elseif(c_type .eq. 'ta')then
            var_2d = 'TAD'
            ext = 'lsx'
            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100,i4time_p
     1w,
     1              ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                     ,field_2d,0,istatus)

            IF(istatus .ne. 1)THEN
                write(6,*)' Error Reading Surface ',var_2d
                goto1200
            endif

            c33_label = 'LAPS Sfc Temp Adv (x 1e-5 Dg K/s)'

            clow = -100.
            chigh = +100.
            cint = 5.

            call make_fnam_lp(i4time_pw,asc9_tim_t,istatus)

            call plot_cont(field_2d,1e-5,clow,chigh,cint,
     1  asc9_tim_t,c33_label,i_overlay,c_display,'nest7grid'
     1                                          ,lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

        elseif(c_type .eq. 'th')then
            var_2d = 'TH'
            ext = 'lsx'
            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100,i4time_p
     1w,
     1              ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                     ,field_2d,0,istatus)

            IF(istatus .ne. 1)THEN
                write(6,*)' Error Reading Surface ',var_2d
                goto1200
            endif

            c33_label = 'LAPS Sfc Potential Temp   (Deg K)'

            clow = +240.
            chigh = +320.
            cint = 2.

            call make_fnam_lp(i4time_pw,asc9_tim_t,istatus)

            call plot_cont(field_2d,1e-0,clow,chigh,cint,
     1  asc9_tim_t,c33_label,i_overlay,c_display,'nest7grid'
     1                                          ,lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

        elseif(c_type .eq. 'te')then
            var_2d = 'THE'
            ext = 'lsx'
            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100,i4time_p
     1w,
     1              ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                     ,field_2d,0,istatus)

            IF(istatus .ne. 1)THEN
                write(6,*)' Error Reading Surface ',var_2d
                goto1200
            endif

            c33_label = 'LAPS Sfc Equiv Potl Temp  (Deg K)'

            clow = +240.
            chigh = +350.
            cint = 2.

            call make_fnam_lp(i4time_pw,asc9_tim_t,istatus)

            call plot_cont(field_2d,1e-0,clow,chigh,cint,
     1  asc9_tim_t,c33_label,i_overlay,c_display,'nest7grid'
     1                                          ,lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

        elseif(c_type .eq. 'vo')then
            var_2d = 'VOR'
            ext = 'lsx'
            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100,i4time_p
     1w,
     1              ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                     ,field_2d,0,istatus)

            IF(istatus .ne. 1)THEN
                write(6,*)' Error Reading Surface ',var_2d
                goto1200
            endif

            c33_label = 'LAPS Sfc Vorticity  (x 1e-5 s^-1)'

            clow = -100.
            chigh = +100.
            cint = 5.

            call make_fnam_lp(i4time_pw,asc9_tim_t,istatus)

            call plot_cont(field_2d,1e-5,clow,chigh,cint,
     1  asc9_tim_t,c33_label,i_overlay,c_display,'nest7grid'
     1                                          ,lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

        elseif(c_type .eq. 'mr')then
            var_2d = 'MR'
            ext = 'lsx'
            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100,i4time_p
     1w,
     1              ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                     ,field_2d,0,istatus)

            IF(istatus .ne. 1)THEN
                write(6,*)' Error Reading Surface ',var_2d
                goto1200
            endif

            c33_label = 'LAPS Sfc Mixing Ratio      (g/kg)'

            clow = -100.
            chigh = +100.
            cint = 2.

            call make_fnam_lp(i4time_pw,asc9_tim_t,istatus)

            call plot_cont(field_2d,1e-0,clow,chigh,cint,
     1  asc9_tim_t,c33_label,i_overlay,c_display,'nest7grid'
     1                                          ,lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

        elseif(c_type .eq. 'dv')then
            var_2d = 'DIV'
            ext = 'lsx'
            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100,i4time_p
     1w,
     1              ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                     ,field_2d,0,istatus)

            IF(istatus .ne. 1)THEN
                write(6,*)' Error Reading Surface ',var_2d
                goto1200
            endif

            c33_label = 'LAPS Sfc Divergence  (x 1e-5 s-1)'

            clow = -100.
            chigh = +100.
            cint = 5.

            call make_fnam_lp(i4time_pw,asc9_tim_t,istatus)

            call plot_cont(field_2d,1e-5,clow,chigh,cint,
     1  asc9_tim_t,c33_label,i_overlay,c_display,'nest7grid'
     1                                          ,lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

        elseif(c_type .eq. 'ha')then ! Theta Advection
            var_2d = 'THA'
            ext = 'lsx'
            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100,i4time_p
     1w,
     1              ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                     ,field_2d,0,istatus)

            IF(istatus .ne. 1)THEN
                write(6,*)' Error Reading Surface ',var_2d
                goto1200
            endif

            c33_label = 'LAPS Sfc Theta Adv   (x 1e-5 K/s)'

            clow = -100.
            chigh = +100.
            cint = 5.

            call make_fnam_lp(i4time_pw,asc9_tim_t,istatus)

            call plot_cont(field_2d,1e-5,clow,chigh,cint,
     1  asc9_tim_t,c33_label,i_overlay,c_display,'nest7grid'
     1                                          ,lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

        elseif(c_type .eq. 'ma')then
            var_2d = 'MRA'
            ext = 'lsx'
            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100,i4time_p
     1w,
     1              ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                     ,field_2d,0,istatus)

            IF(istatus .ne. 1)THEN
                write(6,*)' Error Reading Surface ',var_2d
                goto1200
            endif

            c33_label = 'LAPS Sfc Mstr Adv (x 1e-5 g/kg/s)'

            clow = -100.
            chigh = +100.
            cint = 5.

            call make_fnam_lp(i4time_pw,asc9_tim_t,istatus)

            call plot_cont(field_2d,1e-5,clow,chigh,cint,
     1  asc9_tim_t,c33_label,i_overlay,c_display,'nest7grid'
     1                                          ,lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

        elseif(c_type .eq. 'sp')then
            ext = 'lsx'
            var_2d = 'U'
            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100,i4time_p
     1w,
     1              ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                     ,u_2d,0,istatus)
            var_2d = 'V'
            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100,i4time_p
     1w,
     1              ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                     ,v_2d,0,istatus)


            IF(istatus .ne. 1)THEN
                write(6,*)' Error Reading Surface ',var_2d
                goto1200
            endif

            c33_label = 'LAPS Sfc Wind Speed          (kt)'

            clow = 0.
            chigh = +100.
            cint = 10.

            do i = 1,NX_L
            do j = 1,NY_L
                    if(u_2d(i,j) .eq. r_missing_data
     1    .or. v_2d(i,j) .eq. r_missing_data)then
                        dir(i,j)  = r_missing_data
                        spds(i,j) = r_missing_data
                    else
                        call uvgrid_to_disptrue(u_2d(i,j),
     1                                  v_2d(i,j),
     1                                  dir(i,j),
     1                                  spds(i,j),
     1                                  lon(i,j)     )
                        spds(i,j) = spds(i,j) / mspkt
                    endif
            enddo ! j
            enddo ! i

            call make_fnam_lp(i4time_pw,asc9_tim_t,istatus)

            call plot_cont(spds,1.,clow,chigh,cint,
     1  asc9_tim_t,c33_label,i_overlay,c_display,'nest7grid'
     1                                          ,lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

        elseif(c_type .eq. 'cs')then
            var_2d = 'CSS'
            ext = 'lsx'
            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100,i4time_p
     1w,
     1              ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                     ,field_2d,0,istatus)

            IF(istatus .ne. 1)THEN
                write(6,*)' Error Reading Surface ',var_2d
                goto1200
            endif

            c33_label = 'LAPS  Colorado Severe Storm Index'

            clow = 0.
            chigh = +100.
            cint = 10.

            call make_fnam_lp(i4time_pw,asc9_tim_t,istatus)

            call plot_cont(field_2d,1e-0,clow,chigh,cint,
     1  asc9_tim_t,c33_label,i_overlay,c_display,'nest7grid'
     1                                          ,lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

        elseif(c_type .eq. 'vs')then
            var_2d = 'VIS'
            ext = 'lsx'
            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100,i4time_p
     1w,
     1              ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                     ,field_2d,0,istatus)

            IF(istatus .ne. 1)THEN
                write(6,*)' Error Reading Surface ',var_2d
                goto1200
            endif

            c33_label = 'LAPS Sfc Visibility       (miles)'

            clow = 0.
            chigh = +100.
            cint = -0.1

            call make_fnam_lp(i4time_pw,asc9_tim_t,istatus)

            call plot_cont(field_2d,1600.,clow,chigh,cint,
     1  asc9_tim_t,c33_label,i_overlay,c_display,'nest7grid'
     1                                          ,lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

        elseif(c_type .eq. 'fw')then
            var_2d = 'FWX'
            ext = 'lsx'
            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100,i4time_p
     1w,
     1              ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                     ,field_2d,0,istatus)

            IF(istatus .ne. 1)THEN
                write(6,*)' Error Reading Surface ',var_2d
                goto1200
            endif

            c33_label = 'LAPS Fire Weather          (0-20)'

            clow = 0.
            chigh = +20.
            cint = 2.0

            call make_fnam_lp(i4time_pw,asc9_tim_t,istatus)

            call plot_cont(field_2d,1.,clow,chigh,cint,
     1  asc9_tim_t,c33_label,i_overlay,c_display,'nest7grid'
     1                                          ,lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

        elseif(c_type .eq. 'sc')then
            var_2d = 'SC'
            ext = 'lm2'
            call get_laps_2dgrid(i4time_ref,laps_cycle_time*100,i4time_p
     1w,
     1              ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                     ,field_2d,0,istatus)

            IF(istatus .ne. 1)THEN
                write(6,*)' Error Reading Snow Cover'
                goto1200
            endif

!           c33_label = 'LAPS Snow Cover       (PERCENT)  '
            c33_label = 'LAPS Snow Cover       (TENTHS)   '

            clow = 0.
!           chigh = +100.
!           cint = 20.
            chigh = +10.
            cint = 2.

            call make_fnam_lp(i4time_pw,asc9_tim_t,istatus)

!           call plot_cont(field_2d,1e-2,clow,chigh,cint,
            call plot_cont(field_2d,1e-1,clow,chigh,cint,
     1  asc9_tim_t,c33_label,i_overlay,c_display,'nest7grid'
     1                                          ,lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

        elseif(c_type .eq. 'cb' .or. c_type .eq. 'cc')then

            if(c_type .eq. 'cb')then ! Cloud Base
                ext = 'lcb'
                var_2d = 'LCB'
                call get_laps_2dgrid(i4time_ref,
     1              laps_cycle_time*100,i4time_nearest,
     1              ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                  ,cloud_ceil,0,istatus)

                IF(istatus .ne. 1 .and. istatus .ne. -1)THEN
                    write(6,*)' Error Reading Cloud Base'
                    goto1200
                endif

                c33_label = 'LAPS Cloud Base         m   MSL  '
                clow = 0.
                chigh = 10000.
                cint = 1000.

            elseif(c_type .eq. 'cc')then ! Cloud Ceiling

                ext = 'lcb'
                var_2d = 'CCE'
                call get_laps_2dgrid(i4time_ref,
     1              laps_cycle_time*100,i4time_nearest,
     1              ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                  ,cloud_ceil,0,istatus)

                IF(istatus .ne. 1 .and. istatus .ne. -1)THEN
                    write(6,*)' Error Reading Cloud Ceiling'
                    goto1200
                endif

                c33_label = 'LAPS Cloud Ceiling      m   AGL  '
                clow = 0.
                chigh = 0.
                cint = -100.

            endif

            call make_fnam_lp(i4time_nearest,asc9_tim,istatus)

            call plot_cont(cloud_ceil,1e0,
     1               clow,chigh,cint,asc9_tim,c33_label,
     1          i_overlay,c_display,'nest7grid',lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

        elseif(c_type .eq. 'ct')then
            var_2d = 'LCT'
            ext = 'lcb'
            call get_laps_2dgrid(i4time_ref,864000,i4time_nearest,
     1              ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                     ,cloud_top,0,istatus)

            IF(istatus .ne. 1 .and. istatus .ne. -1)THEN
                write(6,*)' Error Reading Cloud Top'
                goto1200
            endif

            call make_fnam_lp(i4time_nearest,asc9_tim,istatus)

            clow = 0.
            chigh = 20000.
            cint = 1000.
            call plot_cont(cloud_top,1e0,
     1               clow,chigh,cint,asc9_tim,c33_label,
     1          i_overlay,c_display,'nest7grid',lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

        elseif(c_type .eq. 'cv')then
            write(6,2514)
2514        format('     Enter Level (1-42); [-bbbb] for mb; '
     1                          ,'OR [0] for max in column',5x,'? '$)
            read(lun,*)k_level

            if(k_level .gt. 0)then
                var_2d = 'lc3'
                ext = 'lc3'

                call get_laps_2dgrid(i4time_ref,86400,i4time_nearest,
     1           ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                  ,cloud_cvr,k_level,istatus)

            elseif(k_level .lt. 0)then ! k_level is -pressure in mb
               var_2d = 'lcp'
               ext = 'lcp'

               call get_laps_2dgrid(i4time_ref,86400,i4time_nearest,
     1           ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                  ,cloud_cvr,-k_level,istatus)

            else ! k_level .eq. 0
               var_2d = 'lcv'
               ext = 'lcv'

               call get_laps_2dgrid(i4time_ref,86400,i4time_nearest,
     1           ext,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                                  ,cloud_cvr,-k_level,istatus)

            endif

            call make_fnam_lp(i4time_nearest,asc9_tim,istatus)

!           Get Cloud Cover
            if(k_level .gt. 0)then
                read(comment_2d,3515)cloud_height
3515            format(e20.8)

                write(c33_label,3516)nint(cloud_height)
3516            format('LAPS ',i5,'  M MSL   Cloud Cover  ')

                write(6,*)' LVL_CLD = ',lvl_cld

            elseif(k_level .eq. 0)then
                c33_label = 'LAPS Cloud Cover                 '

            else ! k_level .lt. 0
                write(c33_label,3517)-k_level
3517            format('LAPS ',i5,'  MB    Cloud Cover    ')

            endif

            clow = 0.2
            chigh = 0.8
            cint = 0.2
            call plot_cont(cloud_cvr,1e0,
     1               clow,chigh,cint,asc9_tim,c33_label,
     1          i_overlay,c_display,'nest7grid',lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

        elseif(c_type .eq. 'tr')then
            clow = -400.
            chigh = +5000.
            cint = +200.
            c33_label = '                                 '
            asc9_tim_t = '         '
            call plot_cont(topo,1e0,
     1               clow,chigh,cint,asc9_tim_t,c33_label,
     1          i_overlay,c_display,'nest7grid',lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

            i4time_topo = 0

        elseif(c_type .eq. 'gr')then
            call plot_grid(i_overlay,c_display,lat,lon,
     1                     NX_L,NY_L,laps_cycle_time)

        elseif(c_type .eq. 'lf')then
            clow = .5
            chigh = .5
            cint = .5
            c33_label = '                                 '
            asc9_tim_t = '         '
            call plot_cont(rlaps_land_frac,1e0,
     1               clow,chigh,cint,asc9_tim_t,c33_label,
     1          i_overlay,c_display,'nest7grid',lat,lon,jdot,
     1  NX_L,NY_L,r_missing_data,laps_cycle_time)

            i4time_topo = 0

        elseif(c_type .eq. 'cf')then
            call frame
            close(8)

        elseif(c_type .eq. 'q ')then
            goto9000

        endif ! c_field

        goto1200

9000    if(c_display .eq. 'm' .or. c_display .eq. 'p')then
            call frame
        else
            if(c_display .eq. 't')then
            elseif(c_display .eq. 'r')then
                call frame2(c_display)
            else
                call frame
            endif
        endif

        return
        end


        subroutine plot_cont(array,scale,clow,chigh,cint,
     1    asc_tim_9,c33_label,i_overlay,c_display,c_file,lat,lon,jdot,
     1    NX_L,NY_L,r_missing_data,laps_cycle_time)

!       97-Aug-14     Ken Dritz     Added NX_L, NY_L as dummy arguments
!       97-Aug-14     Ken Dritz     Added r_missing_data, laps_cycle_time
!                                   as dummy arguments
!       97-Aug-14     Ken Dritz     Changed LAPS_DOMAIN_FILE to 'nest7grid'
!       97-Aug-14     Ken Dritz     Removed include of lapsparms.for

        common /MCOLOR/mini,maxi

        real*4 lat(NX_L,NY_L),lon(NX_L,NY_L)

        character c33_label*33,asc_tim_9*9,c_metacode*2,asc_tim_24*24
        character*1 c_display
        character*(*) c_file

        real*4 array(NX_L,NY_L)
        real*4 array_plot(NX_L,NY_L)

!       integer*4 ity,ily,istatus
!       data ity/35/,ily/1010/

        include 'icolors.inc'

        Y_SPACING = 3

        write(6,1505)c33_label,scale,asc_tim_9
1505    format(7x,a33,4x,'Units = ',1pe9.0,6x,a9)

        if(asc_tim_9 .ne. '         ')then
            call i4time_fname_lp(asc_tim_9,I4time_file,istatus)
            call cv_i4tim_asc_lp(i4time_file,asc_tim_24,istatus)
            asc_tim_24 = asc_tim_24(1:14)//asc_tim_24(16:17)//' '
        else
            asc_tim_24 = '                        '
        endif


        vmax = -1e30
        vmin = 1e30

        do i = 1,NX_L
        do j = 1,NY_L
            if(array(i,j) .ne. r_missing_data)then
                array_plot(i,j) = array(i,j) / scale
            else
                array_plot(i,j) = array(i,j) 
            endif
            vmax = max(vmax,array_plot(i,j))
            vmin = min(vmin,array_plot(i,j))
        enddo ! i
        enddo ! j

        X_SPACING = NX_L / 28

        do j=NY_L,1,-Y_SPACING
            write(6,500)
     1  (nint(min(max(array_plot(i,j),-99.),999.)),i=1,NX_L,X_SPACING)
500         format(1x,42i3)
        enddo ! j

!       Set Map Background stuff
        if(c_display .eq. 'r' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'l' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 't' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'e' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'n' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'o' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'c' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'v')then
            goto990
        elseif(c_display .eq. 'p')then ! Generate a Map background only
            c_metacode = 'm'
            call lapsplot(array_plot,NX_L,NY_L,clow,chigh,cint,lat,lon
     1          ,c_metacode,c_file,'nest7grid',jdot)
            goto990
        else
            c_metacode = 'c'
        endif

        if(c_metacode .eq. 'c ')then
            i_overlay = i_overlay + 1

        else if(c_metacode .eq. 'm ')then
            write(6,*)' c_metacode,i_overlay = ',c_metacode,i_overlay

            if(c_display .eq. 'r')then
                call lapsplot(array_plot,NX_L,NY_L,clow,chigh,cint,lat,l
     1on
     1          ,c_metacode,c_file,'nest7grid',jdot)
            endif

            c_metacode = 'c '
            i_overlay = 1

            i4time_plot = i4time_file/laps_cycle_time*laps_cycle_time
!       1                                            -laps_cycle_time
            call setusv_dum(2hIN,34)

            iflag = 0

            call get_maxstns(maxstns,istatus)
            if (istatus .ne. 1) then
               write (6,*) 'Error getting value of maxstns'
               stop
            endif

            call plot_station_locations(i4time_plot,lat,lon,NX_L,NY_L,if
     1lag,maxstns)
        endif

        write(6,*)' Plotting: c_metacode,i_overlay = ',
     1                          c_metacode,i_overlay

        if(clow .eq. 0. .and. chigh .eq. 0. .and. cint .gt. 0.)then
            clow =  (nint(vmin/cint)-1) * cint
            chigh = (nint(vmax/cint)+1) * cint
        endif

        write(6,*)' CLOW,HIGH,CINT ',clow,chigh,cint
        write(6,*)' Max/Min = ',vmax,vmin

        call setusv_dum(2hIN,icolors(i_overlay))

        if(c_metacode .ne. 'n ')then
            if(c_metacode .eq. 'c ')then
!                c33_label = 'FSL/'//c33_label(1:29)
                 call upcase(c33_label,c33_label)
                 call set(.00,1.0,.00,1.0,.00,1.0,.00,1.0,1)
                 call write_label_lplot(NX_L,NY_L,c33_label,asc_tim_24
     1                                                      ,i_overlay)
            endif

            if(c_display .ne. 't')then
                mini = icolors(i_overlay)
                maxi = icolors(i_overlay)
            else
                mini = icolors_tt(i_overlay)
                maxi = icolors_tt(i_overlay)
            endif

            call lapsplot(array_plot,NX_L,NY_L,clow,chigh,cint,lat,lon
     1          ,c_metacode,c_file,'nest7grid',jdot)
        endif

990     return
        end



        subroutine plot_barbs(u,v,lat,lon,topo,size,interval,asc_tim_9,
     1  c33_label,
     1  c_field,k_level,i_overlay,c_display,c_file,imax,jmax,kmax,
     1  grid_ra_ref,grid_ra_vel,NX_L,NY_L,r_missing_data,
     1  laps_cycle_time,jdot)      

!       97-Aug-14     Ken Dritz     Added NX_L, NY_L as dummy arguments
!       97-Aug-14     Ken Dritz     Added r_missing_data, laps_cycle_time
!                                   as dummy arguments
!       97-Aug-14     Ken Dritz     Changed LAPS_DOMAIN_FILE to 'nest7grid'
!       97-Aug-14     Ken Dritz     Removed include of lapsparms.for

        character c33_label*33,asc_tim_9*9,c_metacode*2,asc_tim_24*24
        character c_field*2,c_display*1
        character*(*) c_file

        real*4 u(NX_L,NY_L)
        real*4 v(NX_L,NY_L)
        real*4 lat(NX_L,NY_L)
        real*4 lon(NX_L,NY_L)
        real*4 topo(NX_L,NY_L)

        real*4 grid_ra_ref(imax,jmax,kmax)
        real*4 grid_ra_vel(imax,jmax,kmax)

!       integer*4 ity,ily,istatus
!       data ity/35/,ily/1010/

        include 'icolors.inc'

        logical l_obs

        write(6,1505)c33_label,asc_tim_9
1505    format(2x,a33,2x,a9)

        call i4time_fname_lp(asc_tim_9,I4time_file,istatus)
        call cv_i4tim_asc_lp(i4time_file,asc_tim_24,istatus)
        asc_tim_24 = asc_tim_24(1:14)//asc_tim_24(16:17)//' '

!       Set Map Background stuff
        if(c_display .eq. 'r' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'l' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 't' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'e' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'n' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'o' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'c' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'v')then
            goto990
        else
            c_metacode = 'c'
        endif

        if(c_metacode .eq. 'c ')then
            i_overlay = i_overlay + 1
        else if(c_metacode .eq. 'm ')then
            write(6,*)' c_metacode,i_overlay = ',c_metacode,i_overlay

            if(c_display .eq. 'r')then
                call lapsplot(array_plot,NX_L,NY_L,clow,chigh,cint,lat,l
     1on
     1          ,c_metacode,c_file,'nest7grid',jdot)
            endif

            c_metacode = 'c '
            i_overlay = 1

            i4time_plot = i4time_file/laps_cycle_time*laps_cycle_time
!       1                                            -laps_cycle_time
            call setusv_dum(2hIN,34)

            iflag = 0

            call get_maxstns(maxstns,istatus)
            if (istatus .ne. 1) then
               write (6,*) 'Error getting value of maxstns'
               stop
            endif

            call plot_station_locations(i4time_plot,lat,lon,NX_L,NY_L,if
     1lag,maxstns)
        endif

        write(6,*)' Plotting: c_metacode,i_overlay = ',
     1                          c_metacode,i_overlay

        call setusv_dum(2hIN,icolors(i_overlay))

        if(c_metacode .ne. 'n ')then
            if(c_metacode .eq. 'y '
     1  .or. c_metacode .eq. 'c ')then
                 call upcase(c33_label,c33_label)
                 call set(.00,1.0,.00,1.0,.00,1.0,.00,1.0,1)
                 call write_label_lplot(NX_L,NY_L,c33_label,asc_tim_24
     1                                                      ,i_overlay)
            endif


            if(c_metacode .eq. 'y ' .or. c_metacode .eq. 'c ')then
                if(c_field .eq. 'ob')then
                    call getset(mxa,mxb,mya,myb,umin,umax,vmin,vmax,ltyp
     1e)
                    write(6,2031)
2031                format('         Enter Radar #   ',45x,'? '$)
                    read(5,*)i_radar

                    call plot_obs(k_level,.true.,asc_tim_9(1:7)//'00'
     1                  ,i_radar,imax,jmax,kmax
     1                  ,grid_ra_ref,grid_ra_vel,lat,lon,topo,1)
                    return
                endif

                call setusv_dum(2hIN,icolors(i_overlay))

                call get_border(NX_L,NY_L,x_1,x_2,y_1,y_2)

                call set(x_1,x_2,y_1,y_2,1.,float(NX_L),1.,float(NY_L))

                call plot_winds_2d(u,v,interval,size
     1          ,NX_L,NY_L,lat,lon,r_missing_data)
!               call frame

            endif

        endif

990     return
        end


        subroutine plot_grid(i_overlay,c_display,lat,lon,
     1                       NX_L,NY_L,laps_cycle_time)

!       97-Aug-14     Ken Dritz     Added NX_L, NY_L as dummy arguments
!       97-Aug-14     Ken Dritz     Added laps_cycle_time as dummy argument
!       97-Aug-14     Ken Dritz     Changed LAPS_DOMAIN_FILE to 'nest7grid'
!       97-Aug-14     Ken Dritz     Removed include of lapsparms.for

        character c33_label*33,asc_tim_9*9,c_metacode*2,asc_tim_24*24
        character c_field*2,c_display*1

        real*4 lat(NX_L,NY_L)
        real*4 lon(NX_L,NY_L)

!       integer*4 ity,ily,istatus
!       data ity/35/,ily/1010/

        include 'icolors.inc'

        logical l_obs

!       Set Map Background stuff
        if(c_display .eq. 'r' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'l' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 't' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'e' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'n' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'o' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'c' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'v')then
            goto990
        else
            c_metacode = 'c'
        endif

        if(c_metacode .eq. 'c ')then
            i_overlay = i_overlay + 1
        else if(c_metacode .eq. 'm ')then
            write(6,*)' c_metacode,i_overlay = ',c_metacode,i_overlay

            if(c_display .eq. 'r')then
                call lapsplot(array_plot,NX_L,NY_L,clow,chigh,cint
     1                       ,lat,lon
     1                       ,c_metacode,c_file,'nest7grid',jdot)
            endif

            c_metacode = 'c '
            i_overlay = 1

            i4time_plot = i4time_file/laps_cycle_time*laps_cycle_time
!       1                                            -laps_cycle_time
            call setusv_dum(2hIN,34)

        endif

        write(6,*)' Plotting: c_metacode,i_overlay = ',
     1                          c_metacode,i_overlay

        call setusv_dum(2hIN,icolors(i_overlay))

        if(c_metacode .ne. 'n ')then
            if(c_metacode .eq. 'y '
     1  .or. c_metacode .eq. 'c ')then
                 call upcase(c33_label,c33_label)
                 call set(.00,1.0,.00,1.0,.00,1.0,.00,1.0,1)
!                call pwrity(cpux(320),cpux(ity),c33_label,33,2,0,0)      
!                call pwrity
!    1                (cpux(800),cpux(ity),asc_tim_24(1:17),17,2,0,0)
                 call write_label_lplot(NX_L,NY_L,c33_label,asc_tim_24
     1                                                      ,i_overlay)
            endif


            if(c_metacode .eq. 'y ' .or. c_metacode .eq. 'c ')then

                call setusv_dum(2hIN,icolors(i_overlay))

!               call get_border(NX_L,NY_L,x_1,x_2,y_1,y_2)

!               call set(x_1,x_2,y_1,y_2,1.,float(NX_L),1.,float(NY_L))

                call plot_grid_2d(interval,size,NX_L,NY_L,lat,lon)

            endif

        endif

990     return
        end

        subroutine plot_cldpcp_type(cldpcp_type_2d
     1     ,asc_tim_9,c33_label,c_field,k_level,i_overlay,c_display
     1     ,lat,lon,ifield_2d,c_file
     1     ,NX_L,NY_L,laps_cycle_time,jdot)

!       97-Aug-14     Ken Dritz     Added NX_L, NY_L, laps_cycle_time as
!                                   dummy arguments
!       97-Aug-14     Ken Dritz     Removed include of lapsparms.for

        character c33_label*33,asc_tim_9*9,c_metacode*2,asc_tim_24*24
        character c_field*2,c_display*1
        character*(*) c_file

        character cldpcp_type_2d(NX_L,NY_L)
        real*4 lat(NX_L,NY_L)
        real*4 lon(NX_L,NY_L)
        integer*4 ifield_2d(NX_L,NY_L)

        integer*4 iarg

!       integer*4 ity,ily,istatus
!       data ity/35/,ily/1010/

        include 'icolors.inc'

        logical l_obs

        write(6,1505)c33_label,asc_tim_9
1505    format(2x,a33,2x,a9)

        call i4time_fname_lp(asc_tim_9,I4time_file,istatus)
        call cv_i4tim_asc_lp(i4time_file,asc_tim_24,istatus)
        asc_tim_24 = asc_tim_24(1:14)//asc_tim_24(16:17)//' '

!       Set Map Background stuff
        if(c_display .eq. 'r' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'l' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 't' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'e' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'n' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'o' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'c' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'v')then

            interval = 2

            call plot_types_2d(cldpcp_type_2d,interval,size,c_field,.fal
     1se.
     1                                  ,NX_L,NY_L,lat,lon,ifield_2d)
            goto990
        else
            c_metacode = 'c'
        endif

        if(c_metacode .eq. 'c ')then
            i_overlay = i_overlay + 1
        else if(c_metacode .eq. 'm ')then
            write(6,*)' c_metacode,i_overlay = ',c_metacode,i_overlay

            call lapsplot_setup(NX_L,NY_L,lat,lon,jdot)

            c_metacode = 'c '
            i_overlay = 1

            i4time_plot = i4time_file/laps_cycle_time*laps_cycle_time
!       1                                            -laps_cycle_time
            call setusv_dum(2hIN,34)

            iflag = 0

            call get_maxstns(maxstns,istatus)
            if (istatus .ne. 1) then
               write (6,*) 'Error getting value of maxstns'
               stop
            endif

            call plot_station_locations(i4time_plot,lat,lon,NX_L,NY_L,if
     1lag,maxstns)
        endif

        write(6,*)' Plotting: c_metacode,i_overlay = ',
     1                          c_metacode,i_overlay

        call i4time_fname_lp(asc_tim_9,I4time_file,istatus)
        call cv_i4tim_asc_lp(i4time_file,asc_tim_24,istatus)
        asc_tim_24 = asc_tim_24(1:14)//asc_tim_24(16:17)//' '

        call setusv_dum(2hIN,icolors(i_overlay))

        if(c_metacode .ne. 'n ')then
            if(c_metacode .eq. 'y '
     1    .or. c_metacode .eq. 'c ')then
                 call upcase(c33_label,c33_label)
                 call set(.00,1.0,.00,1.0,.00,1.0,.00,1.0,1)
                 call write_label_lplot(NX_L,NY_L,c33_label,asc_tim_24
     1                                                      ,i_overlay)
            endif


            if(c_metacode .eq. 'y ' .or. c_metacode .eq. 'c ')then
                call setusv_dum(2hIN,icolors(i_overlay))

                if(max(NX_L,NY_L) .gt. 61)then
                    interval = 4
                else
                    interval = 2
                endif

                size = 1.0
                call plot_types_2d(cldpcp_type_2d,interval,size,c_field,
     1                       .true.,NX_L,NY_L,lat,lon,ifield_2d)
!               call frame

            endif

        endif

990     return
        end

        subroutine plot_stations(asc_tim_9,c33_label,c_field,i_overlay
     1   ,c_display,lat,lon,c_file,iflag
     1   ,NX_L,NY_L,laps_cycle_time)

!       97-Aug-14     Ken Dritz     Added NX_L, NY_L, laps_cycle_time as
!                                   dummy arguments
!       97-Aug-14     Ken Dritz     Removed include of lapsparms.for

        character c33_label*33,asc_tim_9*9,c_metacode*2,asc_tim_24*24
        character c_field*2,c_display*1
        character*(*) c_file

        real*4 lat(NX_L,NY_L)
        real*4 lon(NX_L,NY_L)

        integer*4 iarg

!       integer*4 ity,ily,istatus
!       data ity/35/,ily/1010/

        include 'icolors.inc'

        logical l_obs

        write(6,1505)c33_label,asc_tim_9
1505    format(2x,a33,2x,a9)

        call i4time_fname_lp(asc_tim_9,I4time_file,istatus)
        call cv_i4tim_asc_lp(i4time_file,asc_tim_24,istatus)
        asc_tim_24 = asc_tim_24(1:14)//asc_tim_24(16:17)//' '

!       Set Map Background stuff
        if(c_display .eq. 'r' .and. i_overlay .eq. 0)then
            c_metacode = 'm'
        elseif(c_display .eq. 'v')then
            goto990
        else
            c_metacode = 'c'
        endif

        if(c_metacode .eq. 'c ')then
            i_overlay = i_overlay + 1
        else if(c_metacode .eq. 'm ')then

            if(iflag .eq. 1)then
                jdot = 0 ! Solid boundaries
            else
                jdot = 1 ! Dotted boundaries
            endif

            write(6,*)' c_metacode,i_overlay = ',c_metacode,i_overlay
            write(6,*)' iflag,jdot = ',iflag,jdot

            call lapsplot_setup(NX_L,NY_L,lat,lon,jdot)

            c_metacode = 'c '
            i_overlay = 1

            i4time_plot = i4time_file ! /laps_cycle_time*laps_cycle_time
!       1                                            -laps_cycle_time
!           call setusv_dum(2hIN,34)
            call setusv_dum(2HIN,11)

            write(6,*)' Call plot_station_locations ',c_metacode

            call get_maxstns(maxstns,istatus)
            if (istatus .ne. 1) then
               write (6,*) 'Error getting value of maxstns'
               stop
            endif

            call plot_station_locations(i4time_plot,lat,lon,NX_L,NY_L,if
     1lag,maxstns)
        endif

        write(6,*)' Plotting: c_metacode,i_overlay = ',
     1                          c_metacode,i_overlay

        call i4time_fname_lp(asc_tim_9,I4time_file,istatus)
        call cv_i4tim_asc_lp(i4time_file,asc_tim_24,istatus)
        asc_tim_24 = asc_tim_24(1:14)//asc_tim_24(16:17)//' '

        call setusv_dum(2hIN,icolors(i_overlay))

        if(c_metacode .ne. 'n ')then
            if(c_metacode .eq. 'y '
     1  .or. c_metacode .eq. 'c ')then
                 call upcase(c33_label,c33_label)
                 call set(.00,1.0,.00,1.0,.00,1.0,.00,1.0,1)
!                call pwrity
!    1                (cpux(320),cpux(ity),c33_label,33,2,0,0)
!                call pwrity
!    1                (cpux(800),cpux(ity),asc_tim_24(1:17),17,2,0,0)
                 call write_label_lplot(NX_L,NY_L,c33_label,asc_tim_24
     1                                                      ,i_overlay)
            endif

            if(c_metacode .eq. 'y ' .or. c_metacode .eq. 'c ')then
                call setusv_dum(2hIN,icolors(i_overlay))

                if(max(NX_L,NY_L) .gt. 61)then
                    interval = 4
                else
                    interval = 2
                endif

                size = 1.0

                write(6,*)' Call plot_station_locations ',c_metacode
                call plot_station_locations(i4time_file,lat,lon
     1                    ,NX_L,NY_L,iflag,maxstns)
            endif

        endif

990     return
        end


        subroutine mklabel33(k_level,c19_label,c33_label)

!       97-Aug-17     Ken Dritz     Lines commented to (temporarily) hardwire
!                                   VERTICAL_GRID at 'PRESSURE' (without
!                                   accessing VERTICAL_GRID)
!       97-Aug-17     Ken Dritz     Removed include of lapsparms.for

        character c19_label*19,c33_label*33

        if(k_level .gt. 0)then
!            if(VERTICAL_GRID .eq. 'HEIGHT')then
!                write(c33_label,101)k_level,c19_label
!101             format('LAPS',I5,' km ',a19)

!            elseif(VERTICAL_GRID .eq. 'PRESSURE')then
                write(c33_label,102)
     1          nint(zcoord_of_level(K_Level)/100.),c19_label
102             format('LAPS',I5,' hPa',a19)

!            endif
        else if(k_level .eq. 0)then
            write(c33_label,103)c19_label
103         format('LAPS  Surface',a19)

        else if(k_level .eq. -1)then
            write(c33_label,104)
104         format('LAPS Steering Winds              ')

        endif

        return
        end

        subroutine plot_station_locations(i4time,lat,lon,ni,nj,iflag,
     1                                    maxstns)

!       97-Aug-14     Ken Dritz     Added maxstns as dummy argument
!       97-Aug-14     Ken Dritz     Removed include of lapsparms.for
!       97-Aug-25     Steve Albers  Removed /read_sfc_cmn/.

!       This routine labels station locations on the H-sect

        real*4 lat(ni,nj),lon(ni,nj)

        real*4 lat_s(maxstns), lon_s(maxstns), elev_s(maxstns)
        real*4 cover_s(maxstns), hgt_ceil(maxstns), hgt_low(maxstns)
        real*4 t_s(maxstns), td_s(maxstns), pr_s(maxstns), sr_s(maxstns)
        real*4 dd_s(maxstns), ff_s(maxstns), ddg_s(maxstns), ffg_s(maxst
     1ns)
        real*4 vis_s(maxstns)
        character stations(maxstns)*3, wx_s(maxstns)*8      ! c5_stamus

!       common /read_sfc_cmn/ lat_s,lon_s,elev_s,cover_s,hgt_ceil,hgt_lo
!    1w
!    1                ,t_s,td_s,pr_s,sr_s,dd_s,ff_s,ddg_s,ffg_s,vis_s
c
        character atime*24, infile*70
        character directory*150,ext*31

        character*9 c9_string
        character*13 filename13

!       Declarations for new read_surface routine
!       New arrays for reading in the SAO data from the LSO files
        real*4   pstn(maxstns),pmsl(maxstns),alt(maxstns),store_hgt(maxs
     1tns,5)
        real*4   ceil(maxstns),lowcld(maxstns),cover_a(maxstns),vis(maxs
     1tns)
     1                                          ,rad(maxstns)

        Integer*4   obstime(maxstns),kloud(maxstns),idp3(maxstns)

        Character   obstype(maxstns)*8
     1             ,store_emv(maxstns,5)*1,store_amt(maxstns,5)*4

        write(6,*)' Reading Station locations from read_sfc for labellin
     1g '
        ext = 'lso'
        call get_directory(ext,directory,len_dir) ! Returns top level directory
        infile = directory(1:len_dir)//filename13(i4time,ext(1:3))

        write(6,*)infile
        call read_surface_old(infile,maxstns,atime,n_meso_g,n_meso_pos,
     &           n_sao_g,n_sao_pos_g,n_sao_b,n_sao_pos_b,n_obs_g,n_obs_p
     1os_g,
     &           n_obs_b,n_obs_pos_b,stations,obstype,lat_s,lon_s,
     &           elev_s,wx_s,t_s,td_s,dd_s,ff_s,ddg_s,
     &           ffg_s,pstn,pmsl,alt,kloud,ceil,lowcld,cover_a,rad,idp3,
     1store_emv,
     &           store_amt,store_hgt,vis,obstime,istatus)

100     write(6,*)'     n_obs_b',n_obs_b

        if(n_obs_b .gt. maxstns .or. istatus .ne. 1)then
            write(6,*)' Too many stations, or no file present'
            istatus = 0
            return
        endif

        size = 0.5
        call getset(mxa,mxb,mya,myb,umin,umax,vmin,vmax,ltype)
        du= float(ni) / 300.

        call get_border(ni,nj,x_1,x_2,y_1,y_2)
        call set(x_1,x_2,y_1,y_2,1.,float(ni),1.,float(nj))
!       call setusv_dum(2HIN,11)

        write(6,*)' Plotting Station locations'
!       Plot Stations
        do i = 1,n_obs_b ! num_sfc
            call latlon_to_rlapsgrid(lat_s(i),lon_s(i),lat,lon
     1                          ,ni,nj,xsta,ysta,istatus)

            if(xsta .lt. 1. .or. xsta .gt. float(ni) .OR.
     1         ysta .lt. 1. .or. ysta .gt. float(nj)          )then       
                    goto80
            endif
!           call supcon(lat_s(i),lon_s(i),usta,vsta)

!           IFLAG = 0        --        Station locations only
!           IFLAG = 1        --        FSL Mesonet only (for WWW)
!           IFLAG = 2        --        All Sfc Obs

            if(iflag .ge. 1)then

!             if(obstype(i) .eq. 'MESO' .or. iflag .eq. 2)then
              if(.true.)then

                if(iflag .eq. 1)call setusv_dum(2HIN,14)

                call pwrity (xsta, ysta-du*3.5, stations(i)(1:3)
     1                    , 3, -1, 0, 0)

                relsize = 1.1

                if(iflag .eq. 1)call setusv_dum(2HIN,11)

                call plot_mesoob(dd_s(i),ff_s(i),ffg_s(i)
     1                 ,t_s(i),td_s(i)
     1                 ,pstn(i),xsta,ysta
     1                 ,lat,lon,ni,nj,relsize,11,iflag)


                if(iflag .eq. 1)call setusv_dum(2HIN,33)

                call line(xsta,ysta+du*0.5,xsta,ysta-du*0.5)
                call line(xsta+du*0.5,ysta,xsta-du*0.5,ysta)


              endif

            else ! Write station location only

                call line(xsta,ysta+du,xsta,ysta-du)
                call line(xsta+du,ysta,xsta-du,ysta)


            endif

80      enddo ! i

        if(iflag .eq. 1)then ! special mesonet label 
            call cv_i4tim_asc_lp(i4time_file,atime,istatus)
            atime = atime(1:14)//atime(16:17)//' '
            ix = 512
            iy = 512
            call pwrity(cpux(ix),cpux(iy),atime(1:17),17,1,0,0)
        endif

        return
        end

        subroutine get_border(ni,nj,x_1,x_2,y_1,y_2)

        if(ni .eq. nj)then
            x_1 = .05
            x_2 = .95
            y_1 = .05
            y_2 = .95
        elseif(ni .lt. nj)then
            ratio = float(ni-1) / float(nj-1)
            x_1 = .50 - .45 * ratio
            x_2 = .50 + .45 * ratio
            y_1 = .05
            y_2 = .95
        elseif(ni .gt. nj)then
            ratio = float(nj-1) / float(ni-1)
            x_1 = .05
            x_2 = .95
            y_1 = .50 - .45 * ratio
            y_2 = .50 + .45 * ratio
        endif

        return
        end



        subroutine setusv_dum(c2_dum,icol_in)

        character*2 c2_dum

        common /icol_index/ icol_current

!       icol = min(icol_in,35)
        icol = icol_in

        write(6,*)' Color # ',icol,icol_in

        call GSTXCI(icol)            
        call GSPLCI(icol)          
        call GSPMCI(icol)           
        call GSFACI(icol)                 

        icol_current = icol

        return
        end


        subroutine write_label_lplot(ni,nj,c33_label,asc_tim_24
     1                                                      ,i_overlay)        

        character*33 c33_label
        character*24 asc_tim_24

        call get_border(ni,nj,x_1,x_2,y_1,y_2)

!       Top label
        y_2 = y_2 + .025

!       Bottom label
        y_1 = y_1 - .025 - .035 * float(i_overlay-1)

        ix = 115
        iy = y_2 * 1024
        call pwrity(cpux(ix),cpux(iy),'NOAA/FSL',8,2,0,0)

        ix = 320
        iy = y_1 * 1024
        call pwrity(cpux(ix),cpux(iy),c33_label,33,2,0,0)

        ix = 800
        iy = y_1 * 1024
        call pwrity(cpux(ix),cpux(iy),asc_tim_24(1:17),17,2,0,0)


        return
        end
