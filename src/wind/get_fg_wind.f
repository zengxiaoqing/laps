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

        subroutine get_fg_wind(
     1          i4time_lapswind,ilaps_cycle_time               ! Input
     1          ,NX_L,NY_L,NZ_L                                ! Input
     1          ,u_mdl_curr,v_mdl_curr                         ! Local/Output
     1          ,u_mdl_diff,v_mdl_diff                         ! Output
     1          ,u_laps_fg,v_laps_fg                           ! Output
     1          ,istatus                                       ! Output
     1                                  )

!       1997 Jun     Ken Dritz     Removed dummy argument NZ_L_MAX (unused)

        logical l_rt     ! Realtime OR Batch Test Mode?
        data l_rt/.true./

        character*31 ext_fg

        dimension u_mdl_prev(NX_L,NY_L,NZ_L),v_mdl_prev(NX_L,NY_L,NZ_L)
        dimension u_mdl_curr(NX_L,NY_L,NZ_L),v_mdl_curr(NX_L,NY_L,NZ_L)
        dimension u_mdl_diff(NX_L,NY_L,NZ_L),v_mdl_diff(NX_L,NY_L,NZ_L)

        dimension u_laps_fg(NX_L,NY_L,NZ_L),v_laps_fg(NX_L,NY_L,NZ_L)

        character*3 var_2d

!       Control which LAPS background can be read in
        if(l_rt)then
            ext_fg = 'ram'  ! Use RAMS first guess
!           ext_fg = 'lw4'  ! Use MAPS first guess (dummy extension)
!           ext_fg = 'lw3'  ! Use LAPS persistance first guess
            itime_start = 0
            itime_stop  = 4
            write(6,*)' 3D option ',ext_fg(1:3)
        else ! Non-Realtime
!           ext_fg = 'lw3'            ! Previous LAPS Wind analysis
!           itime_start = 0
!           itime_stop  = 2
!           write(6,*)' 4D option ',ext_fg(1:3)

            ext_fg = 'lw4'  ! Use MAPS first guess (dummy extension)
            itime_start = 0
            itime_stop  = 4
            write(6,*)' 3D option ',ext_fg(1:3)

        endif

        write(6,*)
        write(6,*)' Obtain wind tendency from model first guesses:'

!  ***  Read in MODEL data   ********************************************
        if(.true.)then

500         write(6,*)
            write(6,*)' Getting Wind Tendency from MODEL:'
            write(6,*)' LAPS EXT: ',ext_fg(1:3)
     1          ,' MDL forecast cycles are ',itime_start,itime_stop


            write(6,*)' Reading u_mdl_prev'
            var_2d = 'U3'
            call get_modelfg_3d(i4time_lapswind-ilaps_cycle_time      
     1          ,var_2d,NX_L,NY_L,NZ_L
     1          ,u_mdl_prev,istatus)
            if(istatus .ne. 1)then
                write(6,*)' Aborting from LAPS Wind Anal'
     1                   ,' - Error reading MODEL Wind'
                return
            endif
            write(6,*)'u_mdl_prev(NX_L/2+1,NY_L/2+1,1) = '
     1                ,u_mdl_prev(NX_L/2+1,NY_L/2+1,1)


            write(6,*)' Reading v_mdl_prev'
            var_2d = 'V3'
            call get_modelfg_3d(i4time_lapswind-ilaps_cycle_time
     1          ,var_2d,NX_L,NY_L,NZ_L
     1          ,v_mdl_prev,istatus)
            if(istatus .ne. 1)then
                write(6,*)' Aborting from LAPS Wind Anal'
     1                   ,' - Error reading MODEL Wind'
                return
            endif
            write(6,*)'v_mdl_prev(NX_L/2+1,NY_L/2+1,1) = '
     1                ,v_mdl_prev(NX_L/2+1,NY_L/2+1,1)


            write(6,*)' Reading u_mdl_curr'
            var_2d = 'U3'
            call get_modelfg_3d(i4time_lapswind,var_2d,NX_L,NY_L,NZ_L
     1          ,u_mdl_curr,istatus)
            if(istatus .ne. 1)then
                write(6,*)' Aborting from LAPS Wind Anal'
     1                   ,' - Error reading MODEL Wind'
                return
            endif
            write(6,*)'u_mdl_curr(NX_L/2+1,NY_L/2+1,1) = '
     1                ,u_mdl_curr(NX_L/2+1,NY_L/2+1,1)


            write(6,*)' Reading v_mdl_curr'
            var_2d = 'V3'
            call get_modelfg_3d(i4time_lapswind,var_2d,NX_L,NY_L,NZ_L
     1               ,v_mdl_curr,istatus)
            if(istatus .ne. 1)then
                write(6,*)' Aborting from LAPS Wind Anal'
     1                   ,' - Error reading MODEL Wind'
                return
            endif
            write(6,*)'v_mdl_curr(NX_L/2+1,NY_L/2+1,1) = '
     1                ,v_mdl_curr(NX_L/2+1,NY_L/2+1,1)


!           Subtract wind field to get time tendency
600         write(6,*)' Subtracting Model Winds to get time tendency'
            do k = 1,NZ_L
            do j = 1,NY_L
            do i = 1,NX_L
                u_mdl_diff(i,j,k) = u_mdl_curr(i,j,k)-u_mdl_prev(i,j,k)
                v_mdl_diff(i,j,k) = v_mdl_curr(i,j,k)-v_mdl_prev(i,j,k)
            enddo ! j
            enddo ! i
            enddo ! k

            write(6,*)'u_mdl_diff(NX_L/2+1,NY_L/2+1,1) = '
     1                ,u_mdl_diff(NX_L/2+1,NY_L/2+1,1)
            write(6,*)'v_mdl_diff(NX_L/2+1,NY_L/2+1,1) = '
     1                ,v_mdl_diff(NX_L/2+1,NY_L/2+1,1)

            I4_elapsed = ishow_timer()

        endif

        write(6,*)
        write(6,*)' Obtain wind from model first guess:'

!       Get previous LAPS analysis for first guess
        if(ext_fg(1:3) .eq. 'ram')then
            istat_persist = 0

        else ! Use persistance LAPS analysis as first guess
            i4time_fg = i4time_lapswind - ilaps_cycle_time

            write(6,*)' Reading in LAPS wind analysis'
     1      ,' from previous cycle for 1st Guess ',ext_fg(1:3)
            write(6,*)' MAPS forecast cycles are ',itime_start
     1                                            ,itime_stop

            call get_uv_3d(i4time_fg,NX_L,NY_L,NZ_L,u_laps_fg,v_laps_fg       
     1                          ,ext_fg,istat_persist)

        endif

        I4_elapsed = ishow_timer()
        if(istat_persist .ne. 1)then
            write(6,*)' No Persistance Winds used: ',ext_fg
            write(6,*)' Using Model Winds for LAPS first guess field'       

            do k = 1,NZ_L
            do j = 1,NY_L
            do i = 1,NX_L
                u_laps_fg(i,j,k) = u_mdl_curr(i,j,k)
                v_laps_fg(i,j,k) = v_mdl_curr(i,j,k)
            enddo ! i
            enddo ! j
            enddo ! k

        else ! Add MAPS time tendency to LAPS first guess
            write(6,*)' Adding MAPS time tendency to LAPS first guess'       

            write(6,*)'u_laps_fg(NX_L/2+1,NY_L/2+1,1) = '
     1                ,u_laps_fg(NX_L/2+1,NY_L/2+1,1)
            write(6,*)'v_laps_fg(NX_L/2+1,NY_L/2+1,1) = '
     1                ,v_laps_fg(NX_L/2+1,NY_L/2+1,1)

            do k = NZ_L,1,-1
            do j = 1,NY_L
            do i = 1,NX_L
                if(u_laps_fg(i,j,k) .ne. 1e-30)then
                    u_laps_fg(i,j,k) = u_laps_fg(i,j,k) 
     1                               + u_mdl_diff(i,j,k)
                    v_laps_fg(i,j,k) = v_laps_fg(i,j,k) 
     1                               + v_mdl_diff(i,j,k)

                else ! Fill in missing data below the terrain
                    u_laps_fg(i,j,k) = u_laps_fg(i,j,k+1)
                    v_laps_fg(i,j,k) = v_laps_fg(i,j,k+1)

                endif
            enddo ! i
            enddo ! j
            enddo ! k

        endif

        write(6,*)'u_laps_fg(NX_L/2+1,NY_L/2+1,1) = '
     1            ,u_laps_fg(NX_L/2+1,NY_L/2+1,1)
        write(6,*)'v_laps_fg(NX_L/2+1,NY_L/2+1,1) = '
     1            ,v_laps_fg(NX_L/2+1,NY_L/2+1,1)

        return
        end


        subroutine get_fg_wind_new(
     1          i4time_lapswind,ilaps_cycle_time               ! Input
     1          ,NX_L,NY_L,NZ_L                                ! Input
     1          ,NTMIN,NTMAX                                   ! Input
     1          ,u_mdl_bkg_4d,v_mdl_bkg_4d                     ! Output
     1          ,u_laps_fg,v_laps_fg                           ! Output
     1          ,istatus                                       ! Output
     1                                  )

!       1998    Steve Albers - FSL

        dimension u_mdl_bkg_4d(NX_L,NY_L,NZ_L,NTMIN:NTMAX)
        dimension v_mdl_bkg_4d(NX_L,NY_L,NZ_L,NTMIN:NTMAX)
        dimension u_laps_fg(NX_L,NY_L,NZ_L),v_laps_fg(NX_L,NY_L,NZ_L)

        character*3 var_2d

        write(6,*)
        write(6,*)' Obtain wind tendency from model first guesses:'

        var_2d = 'U3'

        call get_fg_var(
     1          i4time_lapswind,ilaps_cycle_time               ! Input
     1          ,NX_L,NY_L,NZ_L                                ! Input
     1          ,NTMIN,NTMAX                                   ! Input
     1          ,var_2d                                        ! Input
     1          ,u_mdl_bkg_4d                                  ! Output
     1          ,istatus                                       ! Output
     1                                  )

        if(istatus .ne. 1)then
            write(6,*)' Error processing model first guess for ',var_2d       
            return
        endif

        var_2d = 'V3'

        call get_fg_var(
     1          i4time_lapswind,ilaps_cycle_time               ! Input
     1          ,NX_L,NY_L,NZ_L                                ! Input
     1          ,NTMIN,NTMAX                                   ! Input
     1          ,var_2d                                        ! Input
     1          ,v_mdl_bkg_4d                                  ! Output
     1          ,istatus                                       ! Output
     1                                  )

        if(istatus .ne. 1)then
            write(6,*)' Error processing model first guess for ',var_2d       
            return
        endif

        write(6,*)' Using Model Winds for u/v_laps_fg'       

        call move_3d(u_mdl_bkg_4d(1,1,1,0),u_laps_fg,NX_L,NY_L,NZ_L)       
        call move_3d(v_mdl_bkg_4d(1,1,1,0),v_laps_fg,NX_L,NY_L,NZ_L)

        write(6,*)'u_laps_fg(NX_L/2+1,NY_L/2+1,1) = '
     1            ,u_laps_fg(NX_L/2+1,NY_L/2+1,1)
        write(6,*)'v_laps_fg(NX_L/2+1,NY_L/2+1,1) = '
     1            ,v_laps_fg(NX_L/2+1,NY_L/2+1,1)

        return
        end


        subroutine get_fg_var(
     1          i4time_lapswind,ilaps_cycle_time               ! Input
     1          ,NX_L,NY_L,NZ_L                                ! Input
     1          ,NTMIN,NTMAX                                   ! Input
     1          ,var_2d                                        ! Input
     1          ,var_mdl_bkg_4d                                ! Output
     1          ,istatus                                       ! Output
     1                                  )

!       1999    Steve Albers - FSL

        dimension var_mdl_bkg_4d(NX_L,NY_L,NZ_L,NTMIN:NTMAX)

        character*3 var_2d

        write(6,*)
        write(6,*)' Obtain field tendency from model first guesses:'

!  ***  Read in MODEL data   ********************************************
        do NT = NTMIN,NTMAX

          if(NT .ne. 1)then

            write(6,*)' Reading var_mdl_bkg_4d, var/NT= ',var_2d,NT
            call get_modelfg_3d(i4time_lapswind+ilaps_cycle_time*NT      
     1          ,var_2d,NX_L,NY_L,NZ_L
     1          ,var_mdl_bkg_4d(1,1,1,NT),istatus)
            if(istatus .ne. 1)then
                write(6,*)' Aborting from get_fg_var'
     1                   ,' - Error reading MODEL Data ',var_2d,NT
                return
            endif
            write(6,*)'var_mdl_bkg_4d(NX_L/2+1,NY_L/2+1,1,NT) = '
     1                ,var_mdl_bkg_4d(NX_L/2+1,NY_L/2+1,1,NT)

          else ! NT = 1

!           Note that we can eventually read in the data for NT=1 directly
!           after all obs data is processed "the right way" WRT time
!           The current strategy allows for "NULL" testing of doing improved
!           time interpolation of obs.
            do k = 1,NZ_L
            do j = 1,NY_L
            do i = 1,NX_L
                var_mdl_bkg_4d(i,j,k,1) = 2. * var_mdl_bkg_4d(i,j,k,0) 
     1                                -        var_mdl_bkg_4d(i,j,k,-1)        
            enddo ! i
            enddo ! j
            enddo ! k

          endif ! NT .ne. 1

        enddo ! NT

        I4_elapsed = ishow_timer()

        return
        end

