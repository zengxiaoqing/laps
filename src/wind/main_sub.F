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


        subroutine lapswind_anal  (i4time_sys,
     1              NTMIN,NTMAX,                                     ! I
     1              N_MESO,N_SAO,                                    ! I
     1              heights_3d,                                      ! I
     1              istat_wind)                                      ! O

!       1990         Steve Albers  (Original Version)
!       1992 Apr 30  Steve Albers  (Mod for multi-radar Doppler velocities)
!       1993 Feb 20  Steve Albers  Subroutinize the profiler processing
!       1993 Apr     Steve Albers  Accept RAMS background
!       1993 Nov     Steve Albers  Partial OE interface, 1st guess subroutine
!                                  Eliminate arrays.f from other subroutines
!                                  Eliminate common blocks
!       1994 Oct     Steve Albers  Add i4time_radar_a array
!                                  Move Usfc, Vsfc from LW3 to LWM file
!       1994 Nov 01  Steve Albers  Pass additional arguments into laps_anl
!       1995 Aug     Steve Albers  Some of the U and V arrays are now
!                                  equivalenced to make passing them around
!                                  into a generic Barnes routine easier.
!       1995 Nov     K. Brewster   Changes for nyquist velocity array
!                                  and sounding profiles
!       1995 Nov 28  Steve Albers  iradar_cycle_time = 600 (from 60)
!       1995 Dec  8  Steve Albers  Calls to separate ref and vel radar access
!                                  routines (in getradar.f)
!       1996 Aug     Steve Albers  Mods to allow option to rotate to grid
!                                  north.
!       1997 Feb     Steve Albers  Cleaned up some unneeded arrays/variables
!       1997 Mar     Steve Albers  Removed equiv statements.
!       1997 Jun     Ken Dritz     Removed NZ_L_MAX from the argument list
!                                  in the call to get_fg_wind (it is not
!                                  actually used in get_fg_wind).
!       1997 Jun     Ken Dritz     Changed RETURN statement to STOP statement
!                                  in handling of istatus .ne. 1 after return
!                                  from call to get_domain_laps.
!       1997 Jun     Ken Dritz     Added r_missing_data as dummy argument.
!       1997 Jun     Ken Dritz     Removed parameter declaration of NPTS
!                                  as NX_L * NY_L (because it is unused).
!       1997 Jun     Ken Dritz     Removed parameter declaration of K_3D
!                                  as 2 * NZ_L (because it is unused).
!       1997 Jun     Ken Dritz     Pass N_PIREP in call to rdpirep.
!       1997 Jun     Ken Dritz     Pass NX_L and NY_L a second time to
!                                  vert_wind (once for dummy arguments ni,nj
!                                  and once for newly added dummy arguments
!                                  NX_L_MAX and NY_L_MAX).
!       1997 Jun     Ken Dritz     Pass r_missing_data to vert_wind.
!       1998 Feb     Steve Albers  Call rdpirep with cloud drift winds also

        use mem_namelist ! For these global parameters:
                         ! NX_L,NY_L,nk_laps,max_radars
                         ! ,r_missing_data,i2_missing_data
                         ! ,laps_cycle_time
                         ! 
                         ! Plus a set of wind parameters:

        use mem_grid, ONLY: lat,lon,topo

        use mem_wind

#ifdef SINGLE_EXEC
        use mem_background, only: bgfields%u3, bgfields%v3
#endif

        real, pointer, dimension(:,:,:) :: uanl, vanl, wanl 
        real, pointer, dimension(:,:) :: uanl_sfcitrp,vanl_sfcitrp 

        include 'barnesob.inc'
        type (barnesob) :: obs_point(max_wind_obs)                           

        include 'windparms.inc'

!       Laps Analysis Grids
        real grid_laps_wt(NX_L,NY_L,nk_laps)
        real grid_laps_u(NX_L,NY_L,nk_laps)
        real grid_laps_v(NX_L,NY_L,nk_laps)

        include 'main_sub.inc'

!       Housekeeping
        integer       ss_normal,rtsys_no_data
        parameter(ss_normal        =1, ! success
     1            rtsys_no_data    =3)

        character*3 exts(20),ext_in

        real rlat_radar(max_radars),rlon_radar(max_radars)
     1                               ,rheight_radar(max_radars)
        character*4 radar_name(max_radars)
        real v_nyquist_in(max_radars)
        integer n_vel_grids(max_radars),i4time_radar_a(max_radars)
        Integer       ioffset(max_radars)
        Integer       joffset(max_radars)
        Logical       l_offset_radar 

!       Stuff for call to LAPS analysis
!       real upass1(NX_L,NY_L,nk_laps),vpass1(NX_L,NY_L,nk_laps)

!        real uanl(NX_L,NY_L,nk_laps),vanl(NX_L,NY_L,nk_laps) ! WRT True North
!        real, allocatable, dimension(:,:,:) :: wanl

#ifdef MULTI_EXEC
        real, allocatable, dimension(:,:,:,:) ::
     +                     u_mdl_bkg_4d, v_mdl_bkg_4d
#else
        real, pointer,dimension(:,:,:,:) ::
     +                     u_mdl_bkg_4d, v_mdl_bkg_4d 
#endif

        dimension u_laps_fg(NX_L,NY_L,nk_laps),
     1            v_laps_fg(NX_L,NY_L,nk_laps)

        integer  N_3D_FIELDS
        parameter (N_3D_FIELDS = 3)

!       real outarray_4D(NX_L,NY_L,nk_laps,N_3D_FIELDS)

        character*125 comment_2d,comment_a(2)
        character*10 units_2d,units_a(2)

        character*31 EXT
        character*200 fname

!       character*31 ext_fg

!       Stuff for SFC Winds
!        real uanl_sfcitrp(NX_L,NY_L),vanl_sfcitrp(NX_L,NY_L)

        real, allocatable, dimension(:,:,:,:) :: grid_ra_vel
        real, allocatable, dimension(:,:,:,:) :: grid_ra_nyq
!       real grid_ra_vel(NX_L,NY_L,nk_laps,MAX_RADARS)
!       real grid_ra_nyq(NX_L,NY_L,nk_laps,MAX_RADARS)

        integer idx_radar_a(MAX_RADARS)
        character*31 ext_radar(MAX_RADARS)

        real heights_3d(NX_L,NY_L,nk_laps)
        real heights_1d(nk_laps)

!       Local for laps_anl
        integer n_obs_lvl(nk_laps)

        integer i4_loop_total
        save i4_loop_total
        data i4_loop_total/0/

        istat_wind = 0

        NZ_L = nk_laps
        
        ! Set up pointers to mem_wind arrays so names do not need to be changed
        !    Can not do this with the USE rename list
        uanl => wind%uanl
        vanl => wind%vanl
        uanl_sfcitrp => wind%uanl_sfcitrp
        vanl_sfcitrp => wind%vanl_sfcitrp

        call get_directory('log',fname,len)
        open(15,file=fname(1:len)//'wind_stats.log'
     1      ,status='unknown',err=999)

        write(6,*)' Welcome to the LAPS Wind Analysis'

        if(max_radars .eq. 0 .and. l_use_radial_vel)then
            write(6,*)' WARNING: max_radars = 0, '
     1               ,'setting l_use_radial_vel to .false.'
            l_use_radial_vel = .false.
        endif

!       Get actual grid spacing valid at the gridpoint nearest the center
        icen = NX_L/2 + 1
        jcen = NY_L/2 + 1
        call get_grid_spacing_actual_xy(lat(icen,jcen),lon(icen,jcen)       
     1                        ,grid_spacing_actual_mx
     1                        ,grid_spacing_actual_my
     1                        ,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error return from get_grid_spacing_actual_xy'       
            stop
        endif

        grid_spacing_m = grid_spacing_actual_my

!       write(15,*)' I4time of data = ',i4time_sys
10      continue

        ISTAT = INIT_TIMER()

        n_prods_out = 1 ! (Any number > 0)

!       Housekeeping
        n_prods = 7

        n_msg = 1
        n_pig = 2
        n_prg = 3
        n_sag = 4
        n_d00 = 5
        n_lw3 = 6
        n_lwm = 7

        exts(n_msg) = 'msg'
        exts(n_pig) = 'pig'
        exts(n_prg) = 'prg'
        exts(n_sag) = 'sag'
        exts(n_d00) = 'd00'
        exts(n_lw3) = 'lw3'
        exts(n_lwm) = 'lwm'

        write(6,*)' swlat ',lat(   1,   1)
        write(6,*)' nwlat ',lat(   1,NY_L)
        write(6,*)' selat ',lat(NX_L,   1)
        write(6,*)' nelat ',lat(NX_L,NY_L)

        write(6,*)' swlon ',lon(   1,   1)
        write(6,*)' nwlon ',lon(   1,NY_L)
        write(6,*)' selon ',lon(NX_L,   1)
        write(6,*)' nelon ',lon(NX_L,NY_L)

!       Initialize grids with the missing data value
        do k = 1,NZ_L
        do j = 1,NY_L
        do i = 1,NX_L
            grid_laps_u(i,j,k) = r_missing_data
            grid_laps_v(i,j,k) = r_missing_data
            grid_laps_wt(i,j,k) = r_missing_data
        enddo
        enddo
        enddo

!  ***  Remap Radial Velocities to LAPS Grid  *****************************
        I4_elapsed = ishow_timer()

!       i4_tol = max(laps_cycle_time / 2, iradar_cycle_time / 2)
        i4_tol = 900 ! seconds

        call get_l_offset_radar(nx_l,ny_l,grid_spacing_m,         ! I
     1                          nx_r,ny_r,igrid_r,l_offset_radar) ! O

        allocate(grid_ra_vel(nx_r,ny_r,nz_l,max_radars)
     1           ,STAT=istat_alloc)             
        if(istat_alloc .ne. 0)then
            write(6,*)' ERROR: Could not allocate grid_ra_vel'      
                        stop
        else
            write(6,*)' Allocated grid_ra_vel ',nx_r,ny_r
        endif

        allocate(grid_ra_nyq(nx_r,ny_r,nz_l,max_radars)
     1           ,STAT=istat_alloc)             
        if(istat_alloc .ne. 0)then
            write(6,*)' ERROR: Could not allocate grid_ra_nyq'      
                        stop
        else
            write(6,*)' Allocated grid_ra_nyq ',nx_r,ny_r
        endif

        if(l_use_radial_vel)then
            write(6,*)
            write(6,*)' Reading radial velocity data'
            call get_multiradar_vel(i4time_sys,i4_tol                ! I
     1       ,i4time_radar_a                                         ! O
     1       ,max_radars                                             ! I
     1       ,n_radars                                               ! O
     1       ,ext_radar                                              ! O
     1       ,r_missing_data,NX_L,NY_L,NZ_L,lat,lon                  ! I
     1       ,nx_r,ny_r,igrid_r                                      ! I
     1       ,grid_ra_vel,grid_ra_nyq,idx_radar_a,v_nyquist_in       ! O
     1       ,ioffset,joffset                                        ! O
     1       ,l_offset_radar                                         ! I
     1       ,n_vel_grids                                            ! O
     1       ,rlat_radar,rlon_radar,rheight_radar,radar_name         ! O
     1       ,istat_radar_vel,istat_radar_nyq)                       ! O

            if(n_radars .gt. 0)then
                write(6,*)
                write(6,5545)(radar_name(i),i=1,n_radars)
5545            format(' Retrieved radar Names:',30(1x,a4))
            endif

            if(istat_radar_vel .eq. 1)then
                write(6,*)' Radar 3d vel data successfully read in'
     1                      ,(n_vel_grids(i),i=1,n_radars)
            else
                write(6,*)' Radar 3d vel data NOT successfully read in'
     1                      ,(n_vel_grids(i),i=1,n_radars)
            endif

        else
            n_radars = 0
            istat_radar_vel = 0
            write(6,*)' Not using radar 3d vel, l_use_radial_vel = '
     1               ,l_use_radial_vel

        endif ! l_use_radial_vel

        write(6,*)

!  ***  Access model first guess*********************************************
        I4_elapsed = ishow_timer()

#ifdef MULTI_EXEC
        allocate (  u_mdl_bkg_4d(NX_L,NY_L,NZ_L,NTMIN:NTMAX) )
        allocate (  v_mdl_bkg_4d(NX_L,NY_L,NZ_L,NTMIN:NTMAX) )
         
        call get_fg_wind_new(
     1           i4time_sys,laps_cycle_time                    ! Input
     1          ,NX_L,NY_L,NZ_L                                ! Input
     1          ,NTMIN,NTMAX                                   ! Input
     1          ,u_mdl_bkg_4d,v_mdl_bkg_4d                     ! Output
     1          ,u_laps_fg,v_laps_fg                           ! Output
     1          ,istatus)                                      ! Output 
        if(istatus .ne. 1)then
            write(6,*)' Abort LAPS wind analysis - no first guess info'
            goto 999
        endif
#else
        u_mdl_bkg_4d => bgfields%u3
        v_mdl_bkg_4d => bgfields%v3
        u_laps_fg = bgfields%u3(:,:,:,0)
        v_laps_fg = bgfields%v3(:,:,:,0)
#endif

        write(6,*)' Column of first guess winds'
        write(6,*)'j  v_laps_fg(23,j,13)'
        do j = 1,NY_L
            write(6,5555)j,v_laps_fg(23,j,13)
 5555       format(i4,4f7.2)
        enddo ! j

!  ***  Read in Wind Obs **************************************************

        I4_elapsed = ishow_timer()

!       Initialize observation data structure
        do i = 1,max_wind_obs
            obs_point(i)%l_withhold = .false.
        enddo ! i

        write(6,*)' Calling get_wind_3d_obs'

        nobs_point = 0

        call get_wind_3d_obs(
     1            NX_L,NY_L,NZ_L,                                 ! I
     1            r_missing_data,i2_missing_data,                 ! I
     1            i4time_sys,heights_3d,heights_1d,               ! I
     1            MAX_PR,MAX_PR_LEVELS,weight_prof,l_use_raob,    ! I
     1            l_use_cdw,                                      ! I
     1            N_SAO,N_PIREP,                                  ! I
     1            lat,lon,topo,                                   ! I
     1            NTMIN,NTMAX,                                    ! I
     1            u_mdl_bkg_4d, v_mdl_bkg_4d,                     ! I
     1            grid_laps_u,grid_laps_v,grid_laps_wt,           ! O
     1            max_wind_obs,obs_point,nobs_point,              ! I/O
     1            rlat_radar(1),rlon_radar(1),rheight_radar(1),   ! I
     1            istat_radar_vel,n_vel_grids(1),                 ! I
     1            grid_ra_vel(1,1,1,1),                           ! I
     1            istatus_remap_pro,                              ! O
     1            istatus                )                        ! O

        if(istatus .ne. 1)then
            write(6,*)' Abort LAPS wind analysis'
            goto 999
        endif

        if(istatus_remap_pro .ne. 1)then
            goto 999
        endif

        num_wind_obs = nobs_point

#ifdef MULTI_EXEC
        deallocate (u_mdl_bkg_4d, v_mdl_bkg_4d)
#endif

        I4_elapsed = ishow_timer()

        if((.not. l_grid_north_bkg) .and. l_grid_north_anal)then       
                write(6,*)' Rotating first guess to grid north'

                do k = 1, NZ_L
                do j = 1, NY_L
                do i = 1, NX_L
                    call uvtrue_to_uvgrid(
     1                           u_laps_fg(i,j,k),v_laps_fg(i,j,k)
     1                          ,u_grid   ,v_grid
     1                          ,lon(i,j)           )
                    u_laps_fg(i,j,k) = u_grid
                    v_laps_fg(i,j,k) = v_grid
                enddo
                enddo
                enddo

                I4_elapsed = ishow_timer()

        endif

        if(l_grid_north_anal)then       
                write(6,*)' Rotating obs arrays to grid north'

                do k = 1, NZ_L
                do j = 1, NY_L
                do i = 1, NX_L
                    if(     grid_laps_u(i,j,k) .ne. r_missing_data
     1                .and. grid_laps_v(i,j,k) .ne. r_missing_data )then
                        call uvtrue_to_uvgrid(
     1                           grid_laps_u(i,j,k),grid_laps_v(i,j,k)
     1                          ,u_grid   ,v_grid
     1                          ,lon(i,j)           )
                        grid_laps_u(i,j,k) = u_grid
                        grid_laps_v(i,j,k) = v_grid
                    endif
                enddo
                enddo
                enddo

                write(6,*)' Rotating obs structure to grid north'
                do iob = 1,nobs_point
                    i = obs_point(iob)%i
                    j = obs_point(iob)%j
                    call uvtrue_to_uvgrid(obs_point(iob)%valuef(1)
     1                                   ,obs_point(iob)%valuef(2)
     1                                   ,u_grid   ,v_grid
     1                                   ,lon(i,j)           )
                    obs_point(iob)%valuef(1) = u_grid
                    obs_point(iob)%valuef(2) = v_grid
                enddo ! iob

                I4_elapsed = ishow_timer()

        endif ! l_grid_north

!sms$serial end

        call laps_anl(grid_laps_u,grid_laps_v
     1       ,obs_point,max_wind_obs,nobs_point                         ! I
     1       ,n_radars,istat_radar_vel                                  ! I
     1       ,grid_ra_vel,grid_ra_nyq,v_nyquist_in,idx_radar_a          ! I
!    1       ,upass1,vpass1                                             ! O
     1       ,n_var                                                     ! I
     1       ,uanl,vanl                                                 ! O
     1       ,grid_laps_wt,weight_bkg_const_wind,rms_thresh_wind        ! I/L
     1       ,r0_barnes_max_m,brns_conv_rate_wind                       ! I
     1       ,max_radars                                                ! I
     1       ,n_vel_grids,rlat_radar,rlon_radar,rheight_radar           ! I
     1       ,thresh_2_radarobs_lvl_unfltrd                             ! I
     1       ,thresh_4_radarobs_lvl_unfltrd                             ! I
     1       ,thresh_9_radarobs_lvl_unfltrd                             ! I
     1       ,thresh_25_radarobs_lvl_unfltrd                            ! I
     1       ,u_laps_fg,v_laps_fg                                       ! I/L
     1       ,NX_L,NY_L,NZ_L,lat,lon                                    ! I
     1       ,nx_r,ny_r,ioffset,joffset                                 ! I
     1       ,i4time_sys,grid_spacing_m                                 ! I
     1       ,r_missing_data                                            ! I
     1       ,heights_3d                                                ! I
     1       ,i4_loop_total                                             ! O
     1       ,l_derived_output,l_grid_north_anal,l_3pass                ! I
     1       ,l_correct_unfolding                                       ! I
!    1       ,n_iter_wind
     1       ,weight_cdw,weight_sfc,weight_pirep,weight_prof            ! I
     1       ,weight_radar                                              ! I
     1       ,istatus)                                                  ! O

!sms$serial(default=ignore)  begin              

        if(istatus .ne. 1)then
                write(6,*)' Error in Wind Analysis (laps_anl)'
                goto 999
        endif

        if(l_grid_north .and. .not. l_grid_north_out)then

                I4_elapsed = ishow_timer()

                write(6,*)' Rotating analyzed wind back to true north'
                do k = 1, NZ_L
                do j = 1, NY_L
                do i = 1, NX_L
                        call uvgrid_to_uvtrue(
     1                           uanl(i,j,k),vanl(i,j,k)
     1                          ,u_true   ,v_true
     1                          ,lon(i,j)           )
                    uanl(i,j,k) = u_true
                    vanl(i,j,k) = v_true
                enddo
                enddo
                enddo

                I4_elapsed = ishow_timer()

        endif ! rotate back to true north

        write(6,*)'uanl(NX_L/2+1,NY_L/2+1,1) = '
     1            ,uanl(NX_L/2+1,NY_L/2+1,1)
        write(6,*)'vanl(NX_L/2+1,NY_L/2+1,1) = '
     1            ,vanl(NX_L/2+1,NY_L/2+1,1)

        I4_elapsed = ishow_timer()

#ifdef MULTI_EXEC
!       Need to allocate actual datatype member, can not use wanl pointer for 
!       allocation
        allocate( wind%wanl(NX_L,NY_L,nk_laps), STAT=istat_alloc )
        if(istat_alloc .ne. 0)then
            write(6,*)' ERROR: Could not allocate wind%wanl'
            stop
        endif
#endif
        wanl => wind%wanl

        call wind_post_process(i4time_sys
     1                        ,uanl,vanl                            ! I
     1                        ,wanl                                 ! O
     1                        ,NX_L,NY_L,NZ_L,N_3D_FIELDS           ! I
     1                        ,heights_3d                           ! I
     1                        ,uanl_sfcitrp,vanl_sfcitrp            ! I
     1                        ,topo,lat,lon,grid_spacing_m          ! I
     1                        ,r_missing_data,l_grid_north_out      ! I
     1                        ,istat_lw3)

        I4_elapsed = ishow_timer()

        if(i4_elapsed .gt. 0)then
            pct_loop_total = 
     1          float(i4_loop_total)/float(i4_elapsed) * 100.
        else
            pct_loop_total = 0.
        endif
        write(6,998)           i4_loop_total                 
     1                        ,i4_elapsed - i4_loop_total
     1                        ,i4_elapsed     
     1                        ,pct_loop_total
998     format(' i4time (FINAL) loop/nonloop/total/% = ',3i5,f8.2)

        istat_wind = 1

        write(6,*)' End of LAPS Wind Analysis'

        close(15)
!sms$serial end

 999    continue
   
        deallocate (grid_ra_vel, grid_ra_nyq)

        return

        end

