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
  

        subroutine get_wind_3d_obs(
     1            NX_L,NY_L,NZ_L,                                   ! I
     1            r_missing_data,i2_missing_data,                   ! I
     1            i4time_lapswind,heights_3d,heights_1d,            ! I
     1            MAX_PR,MAX_PR_LEVELS,weight_prof,l_use_raob,      ! I
     1            l_use_cdw,                                        ! I
     1            N_SFC,N_PIREP,                                    ! I
     1            lat,lon,                                          ! I
     1            NTMIN,NTMAX,                                      ! I
     1            u_laps_fg,v_laps_fg,                              ! O
     1            grid_laps_u,grid_laps_v,grid_laps_wt,             ! I/O
     1            max_obs,obs_point,nobs_point,                     ! I/O
     1            rlat_radar,rlon_radar,rheight_radar,              ! I
     1            istat_radar_vel,n_vel_grids,                      ! I
     1            grid_ra_vel,                                      ! I
     1            istatus_remap_pro,                                ! O
     1            istatus                )                          ! O

!       1997 Jun     Ken Dritz       Added NX_L, NY_L, NZ_L as dummy arguments,
!                                    making arrays with those dimensions
!                                    automatic.
!       1997 Jun     Ken Dritz       Added r_missing_data and i2_missing_data
!                                    as dummy arguments.
!       1997 Jun     Ken Dritz       Removed include of 'lapsparms.for'.

        include 'barnesob.inc'
        type (barnesob) obs_point(max_obs)                           

!       LAPS Grid Dimensions
                                                   
        real*4 lat(NX_L,NY_L)                      
        real*4 lon(NX_L,NY_L)                      

!       Sectorized Uniform Wind

!       Profiler Stuff
        real lat_pr(MAX_PR)                        
        real lon_pr(MAX_PR)                        
        character*8 obstype(MAX_PR)
        character*5 c5_name_a(MAX_PR)

!       Profiler Observations

        integer nlevels_obs_pr(MAX_PR)             
        real ob_pr_u (MAX_PR,NZ_L) ! Vertically interpolated Profiler wind
        real ob_pr_v (MAX_PR,NZ_L) ! Vertically interpolated Profiler wind

!       Laps Analysis Grids
        real grid_laps_wt(NX_L,NY_L,NZ_L)                               
        real grid_laps_u(NX_L,NY_L,NZ_L)                                
        real grid_laps_v(NX_L,NY_L,NZ_L)                                

!       Data flags for LAPS analysis
        logical L_profiler
        parameter (L_profiler = .true.)

        dimension u_mdl_bkg_4d(NX_L,NY_L,NZ_L,NTMIN:NTMAX)
        dimension v_mdl_bkg_4d(NX_L,NY_L,NZ_L,NTMIN:NTMAX)

        dimension u_laps_fg(NX_L,NY_L,NZ_L)
        dimension v_laps_fg(NX_L,NY_L,NZ_L)

        real*4 heights_3d(NX_L,NY_L,NZ_L)
        real*4 heights_1d(NZ_L)
        real*4 grid_ra_vel(NX_L,NY_L,NZ_L)

        character*3 ext_in

        logical l_use_raob,l_use_cdw,l_use_all_nontower_lvls

        l_use_all_nontower_lvls = .false.

!  ***  Read in Model First Guess Data  **************************************

        call get_laps_cycle_time(ilaps_cycle_time,istatus)
        if(istatus .eq. 1)then
            write(6,*)' ilaps_cycle_time = ',ilaps_cycle_time
        else
            write(6,*)' Error getting laps_cycle_time'
            return
        endif

        call get_fg_wind_new(
     1          i4time_lapswind,ilaps_cycle_time               ! Input
     1          ,NX_L,NY_L,NZ_L                                ! Input
     1          ,NTMIN,NTMAX                                   ! Input
     1          ,u_mdl_bkg_4d,v_mdl_bkg_4d                     ! Output
     1          ,u_laps_fg,v_laps_fg                           ! Output
     1          ,istatus)                                      ! Output

        if(istatus .ne. 1)then
            write(6,*)' Abort LAPS wind analysis - no first guess info'
            return
        endif

        write(6,*)' Column of first guess winds'
        write(6,*)'j  v_laps_fg(23,j,13)'
        do j = 1,NY_L
            write(6,5555)j,v_laps_fg(23,j,13)
 5555       format(i4,4f7.2)
        enddo ! j


!  ***  Read in Profiler Data  ********************************************

        write(6,*)' Calling read_profiles'

        call read_profiles(
     1            i4time_lapswind,heights_3d,                       ! I
     1            lat_pr,lon_pr,obstype,c5_name_a,                  ! O
     1            lat,lon,                                          ! I
     1            MAX_PR,MAX_PR_LEVELS,                             ! I
     1            l_use_raob,l_use_all_nontower_lvls,               ! I
     1            ob_pr_u , ob_pr_v ,                               ! O
     1            max_obs,obs_point,nobs_point,weight_prof,         ! I/O
     1            nlevels_obs_pr,n_profiles,                        ! O
     1            rlat_radar,rlon_radar,rheight_radar,              ! I
     1            n_vel_grids,                                      ! I
     1            u_mdl_bkg_4d,v_mdl_bkg_4d,NTMIN,NTMAX,            ! I
     1            ilaps_cycle_time,r_missing_data,                  ! I
     1            NX_L,NY_L,NZ_L,                                   ! I
     1            istatus                )                          ! O

        if(istatus .ne. 1)then
            write(6,*)' Abort read_profiles'
            return
        endif

        I4_elapsed = ishow_timer()

! ***   Remapping + Barnes Analysis of Profiler Data in u & v *****************

        call remap_profiles(
     1           ob_pr_u,ob_pr_v                                  ! I
     1          ,grid_laps_u,grid_laps_v,grid_laps_wt             ! I/O
     1          ,max_obs,obs_point,nobs_point                     ! I/O
     1          ,lat,lon                                          ! I
     1          ,NX_L,NY_L,NZ_L,MAX_PR                            ! I
     1          ,nlevels_obs_pr,lat_pr,lon_pr,obstype,n_profiles  ! I
     1          ,c5_name_a                                        ! I
     1          ,r_missing_data                                   ! I
     1          ,weight_prof                                      ! O
     1          ,l_profiler,l_use_all_nontower_lvls               ! I
     1          ,istatus_remap_pro)                               ! O

!  ***  Read in SFC wind data   *******************************************

	call get_maxstns(maxstns,istatus)
	if (istatus .ne. 1) then
	   write (6,*) 'Error obtaining maxstns'
           return
	endif

        call rdsfc(i4time_lapswind,heights_3d                            ! I
     1            ,N_SFC,maxstns                                         ! I
     1            ,lat,lon,r_missing_data                                ! I
     1            ,n_sfc_obs                                             ! O
     1            ,grid_laps_wt,grid_laps_u,grid_laps_v                  ! I/O
     1            ,max_obs,obs_point,nobs_point                          ! I/O
     1            ,NX_L,NY_L,NZ_L                                        ! I
     1            ,istatus)                                              ! O
        if(istatus .ne. 1)then
            write(6,*)
     1     ' Aborting from get_wind_obs - Error reading sfc obs'
            return
!       else
!           i4time_array(n_sag) = i4time_lapswind
!           j_status(n_sag) = ss_normal
        endif

        I4_elapsed = ishow_timer()

!  ***  Read in Pirep data   ********************************************

        n_pirep_obs = 0

        ext_in = 'pin'

        call rdpoint(i4time_lapswind,heights_3d
     1  ,N_PIREP,n_pirep_obs,ext_in
     1  ,NX_L,NY_L,NZ_L
     1  ,u_mdl_bkg_4d,v_mdl_bkg_4d,NTMIN,NTMAX                         ! I
     1  ,lat,lon
!    1  ,pirep_i,pirep_j,pirep_k,pirep_u,pirep_v
     1  ,grid_laps_wt,grid_laps_u,grid_laps_v                          ! I/O
     1  ,max_obs,obs_point,nobs_point                                  ! I/O
     1  ,istatus)                                                      ! O

        if(istatus .ne. 1)then
            write(6,*)
     1       ' Aborting from get_wind_obs - Error reading pireps'
            return
!       else
!           i4time_array(n_pig) = i4time_lapswind
!           j_status(n_pig) = ss_normal
        endif

!  ***  Read in cloud drift wind data   ***********************************

        if(l_use_cdw)then

            ext_in = 'cdw'

            call rdpoint(i4time_lapswind,heights_3d
     1      ,N_PIREP,n_pirep_obs,ext_in
     1      ,NX_L,NY_L,NZ_L
     1      ,u_mdl_bkg_4d,v_mdl_bkg_4d,NTMIN,NTMAX                     ! I
     1      ,lat,lon
!    1      ,pirep_i,pirep_j,pirep_k,pirep_u,pirep_v
     1      ,grid_laps_wt,grid_laps_u,grid_laps_v                      ! I/O
     1      ,max_obs,obs_point,nobs_point                              ! I/O
     1      ,istatus)                                                  ! O

            if(istatus .ne. 1)then
                write(6,*)
     1           ' Aborting from LAPS Wind Anal - Error reading cdw'
                return
!           else
!               i4time_array(n_pig) = i4time_lapswind
!               j_status(n_pig) = ss_normal
            endif

        endif ! l_use_cdw


        return
        end


        subroutine remap_profiles(
     1           ob_pr_u,ob_pr_v                                     ! I
     1          ,grid_laps_u,grid_laps_v,grid_laps_wt                ! I/O
     1          ,max_obs,obs_point,nobs_point                        ! I/O
     1          ,lat,lon                                             ! I
     1          ,ni,nj,nk,MAX_PR                                     ! I
     1          ,nlevels_obs_pr,lat_pr,lon_pr,obstype,n_profiles     ! I
     1          ,c5_name_a
     1          ,r_missing_data                                      ! I
     1          ,weight_prof                                         ! O
     1          ,l_profiler,l_use_all_nontower_lvls                  ! I
     1          ,istatus)                                            ! O

!       Perform horizontal remapping of profile obs onto LAPS grid
!       They have already been vertically interpolated

        include 'barnesob.inc'
        type (barnesob) obs_point(max_obs)                           

!       Profile Observations

        integer nlevels_obs_pr(MAX_PR)
        real*4 lat_pr(MAX_PR)
        real*4 lon_pr(MAX_PR)
        character*8 obstype(MAX_PR)
        character*5 c5_name_a(MAX_PR)
        real*4 ob_pr_u (MAX_PR,nk) ! Vertically interpolated Profile wind
        real*4 ob_pr_v (MAX_PR,nk) ! Vertically interpolated Profile wind

!       Barnes Profile analysis

        real*4 lat(ni,nj),lon(ni,nj)

        real*4 grid_laps_u(ni,nj,nk)
        real*4 grid_laps_v(ni,nj,nk)
        real*4 grid_laps_wt(ni,nj,nk)

        logical l_profiler, l_use_all_nontower_lvls

        write(6,*)
        write(6,*)' Subroutine remap_profiles: # of profiles = '
     1           ,n_profiles

        do i_pr = 1,n_profiles ! MAX_PR

            if(.not. l_use_all_nontower_lvls .and.
     1        obstype(i_pr)(1:5) .ne. 'TOWER'     )then       

                if(nlevels_obs_pr(i_pr) .gt. 0)then
                    call latlon_to_rlapsgrid(lat_pr(i_pr),lon_pr(i_pr)
     1                                      ,lat,lon,ni,nj,ri,rj
     1                                      ,istatus)     
                    if(istatus .ne. 1)then
                        write(6,*)' NOTE... Profile is outside domain'  
     1                           ,i_pr,i_ob,j_ob,' ',obstype(i_pr)

                    else ! inside domain
                        i_ob = nint(ri)
                        j_ob = nint(rj)

                        write(6,*)' Remapping profile',i_pr,i_ob,j_ob
     1                           ,nlevels_obs_pr(i_pr)
     1                           ,' ',obstype(i_pr)       
     1                           ,' ',c5_name_a(i_pr)       

                        do k = 1,nk
                            if(ob_pr_u(i_pr,k) .ne. r_missing_data)then

                                ob_u = ob_pr_u (i_pr,k)
                                ob_v = ob_pr_v (i_pr,k)

!                 ***           Map Observation onto LAPS grid   ***
                                if(l_profiler)then
!                                   grid_laps_u(i_ob,j_ob,k) = ob_u
!                                   grid_laps_v(i_ob,j_ob,k) = ob_v
!                                   grid_laps_wt(i_ob,j_ob,k) = weight_prof
                                    write(6,11)k,ob_u,ob_v
 11                                 format(10x,i4,2f8.1)

!                                   Add to data structure (still is subsampling)
                                    nobs_point = nobs_point + 1
                                    obs_point(nobs_point)%i = i_ob
                                    obs_point(nobs_point)%j = j_ob
                                    obs_point(nobs_point)%k = k
                                    obs_point(nobs_point)%rk = k
                                    obs_point(nobs_point)%valuef(1)=ob_u       
                                    obs_point(nobs_point)%valuef(2)=ob_v
                                    obs_point(nobs_point)%weight = 
     1                                                    weight_prof       
                                    obs_point(nobs_point)%type = 'prof'      

                                endif ! l_profiler

                            endif ! not missing data
                        enddo ! k
                    endif ! istatus (in/out of domain)

                else
                    write(6,*)' Zero levels / outside domain'
     1                       ,i_pr,i_ob,j_ob,' ',obstype(i_pr)

                endif ! data present

            endif ! remapping with vertical interpolation

        enddo ! i_pr

        I4_elapsed = ishow_timer()

        istatus = 1

        return
        end


