
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

        character*125 comment_2d
        character*10 units_2d
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

        logical l_rt     ! Realtime OR Batch Test Mode?
        data l_rt/.true./

        character*31 ext_fg

        dimension u_mdl_bkg_4d(NX_L,NY_L,NZ_L,NTMIN:NTMAX)
        dimension v_mdl_bkg_4d(NX_L,NY_L,NZ_L,NTMIN:NTMAX)
        dimension u_laps_fg(NX_L,NY_L,NZ_L),v_laps_fg(NX_L,NY_L,NZ_L)

        character*125 comment_2d
        character*10 units_2d
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
        do NT = NTMIN,NTMAX

            write(6,*)' Reading u_mdl_bkg_4d ',NT
            var_2d = 'U3'
            call get_modelfg_3d(i4time_lapswind+ilaps_cycle_time*NT      
     1          ,var_2d,NX_L,NY_L,NZ_L
     1          ,u_mdl_bkg_4d(1,1,1,NT),istatus)
            if(istatus .ne. 1)then
                write(6,*)' Aborting from LAPS Wind Anal'
     1                   ,' - Error reading MODEL Wind'
                return
            endif
            write(6,*)'u_mdl_bkg_4d(NX_L/2+1,NY_L/2+1,1,NT) = '
     1                ,u_mdl_bkg_4d(NX_L/2+1,NY_L/2+1,1,NT)


            write(6,*)' Reading v_mdl_bkg_4d ',NT
            var_2d = 'V3'
            call get_modelfg_3d(i4time_lapswind+ilaps_cycle_time*NT      
     1          ,var_2d,NX_L,NY_L,NZ_L
     1          ,v_mdl_bkg_4d(1,1,1,NT),istatus)
            if(istatus .ne. 1)then
                write(6,*)' Aborting from LAPS Wind Anal'
     1                   ,' - Error reading MODEL Wind'
                return
            endif
            write(6,*)'v_mdl_bkg_4d(NX_L/2+1,NY_L/2+1,1,NT) = '
     1                ,v_mdl_bkg_4d(NX_L/2+1,NY_L/2+1,1,NT)

        enddo ! NT

        write(6,*)
        write(6,*)' Obtain wind from model first guess:'

!       Get previous LAPS analysis for first guess
        if(ext_fg(1:3) .eq. 'ram')then
            istat_persist = 0

        else ! Use persistance LAPS analysis as first guess
            i4time_fg = i4time_lapswind - ilaps_cycle_time

            write(6,*)' Reading in LAPS wind analysis'
     1      ,' from previous cycle for 1st Guess ',ext_fg(1:3)
            write(6,*)' Model forecast cycles are ',itime_start
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
                u_laps_fg(i,j,k) = u_mdl_bkg_4d(i,j,k,0)
                v_laps_fg(i,j,k) = v_mdl_bkg_4d(i,j,k,0)
            enddo ! i
            enddo ! j
            enddo ! k

        else ! Add Model time tendency to LAPS first guess
            write(6,*)' Adding Model time tendency to LAPS first guess'      

            write(6,*)'u_laps_fg(NX_L/2+1,NY_L/2+1,1) = '
     1                ,u_laps_fg(NX_L/2+1,NY_L/2+1,1)
            write(6,*)'v_laps_fg(NX_L/2+1,NY_L/2+1,1) = '
     1                ,v_laps_fg(NX_L/2+1,NY_L/2+1,1)

            do k = NZ_L,1,-1
            do j = 1,NY_L
            do i = 1,NX_L
                if(u_laps_fg(i,j,k) .ne. 1e-30)then
                    u_laps_fg(i,j,k) = u_laps_fg(i,j,k) 
     1                               + u_mdl_bkg_4d(i,j,k,0)
     1                               - u_mdl_bkg_4d(i,j,k,-1)
                    v_laps_fg(i,j,k) = v_laps_fg(i,j,k) 
     1                               + v_mdl_bkg_4d(i,j,k,0)
     1                               - v_mdl_bkg_4d(i,j,k,-1)

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

