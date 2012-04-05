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
        subroutine rdpoint(i4time,heights_3d                           ! I
     1  ,N_POINT                                                       ! I
     1  ,n_point_obs                                                   ! I/O
     1  ,ext_in                                                        ! I
     1  ,ni,nj,nk                                                      ! I
     1  ,u_mdl_bkg_4d,v_mdl_bkg_4d,NTMIN,NTMAX                         ! I
     1  ,lat,lon                                                       ! I
!    1  ,point_i,point_j,point_k,point_u,point_v
     1  ,grid_laps_wt,grid_laps_u,grid_laps_v                          ! O
     1  ,max_obs,obs_point,nobs_point                                  ! I/O
     1  ,istatus)                                                      ! O

!      ~1990        Steve Albers  Original Version
!       Modified 2/1993 by Steve Albers to fix check on pirep being in the
!       domain as suggested by Steve Olson of LL.

!       1997 Jun    Ken Dritz     Added N_POINT as dummy argument, making
!                                 arrays dimensioned therewith automatic.
!       1997 Jun    Ken Dritz     Removed include of 'lapsparms.for'.
!       1998        Steve Albers  Called Sequentially for acars winds, 
!                                 then cloud drift winds.
!	2006 Feb    Yuanfu Xie	  Change rk attribute to rk from k_grid
!				  and assign values of ri and rj for
!			       	  obs_point.

!******************************************************************************

        use mem_namelist, ONLY: iwrite_output

        include 'barnesob.inc'
        type (barnesob) :: obs_point(max_obs)                           

!       LAPS Grid Dimensions

        include 'windparms.inc' ! weight_pirep, weight_cdw

        real lat(ni,nj)
        real lon(ni,nj)

!       Point obs

        integer point_i(N_POINT) ! X point coordinates
        integer point_j(N_POINT) ! Y point coordinates
        integer point_k(N_POINT) ! Z point coordinates
        real    point_u(N_POINT) ! u point component
        real    point_v(N_POINT) ! v point component


!       Laps Analysis Grids
        real grid_laps_wt(ni,nj,nk)
        real grid_laps_u(ni,nj,nk)
        real grid_laps_v(ni,nj,nk)

!******************************************************************************

        real heights_3d(ni,nj,nk)
        real pres_3d(ni,nj,nk)

        dimension u_mdl_bkg_4d(ni,nj,nk,NTMIN:NTMAX)
        dimension v_mdl_bkg_4d(ni,nj,nk,NTMIN:NTMAX)

        character*9 asc9_tim_point
        character ext*31, ext_in*3

        logical l_eof, l_geoalt

        real r_missing_data

        n_wisdom = 0

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1) then
          write(6,*) 'Cannot access r_missing_data value'
          write(6,*) 'Aborting read of file: ',ext_in
          return
        endif

        call get_pres_3d(i4time,ni,nj,nk,pres_3d,istatus)
        if(istatus .ne. 1) then
          write(6,*) 'Bad status returned from get_pres_3d'
          write(6,*) 'Aborting read of file: ',ext_in
          return
        endif

!       Open input intermediate data file
        lun_in = 31
        call open_lapsprd_file_read(lun_in,i4time,ext_in,istatus)
        if(istatus .ne. 1)go to 999

        lun_pig = 32
        ext = 'pig'

!       Open output intermediate data file
        if(ext_in .eq. 'pin')then
            if(iwrite_output .ge. 1)then
                call open_lapsprd_file(lun_pig,i4time,ext,istatus)
                if(istatus .ne. 1)go to 888
            endif

            call get_windob_time_window('ACARS',i4_window_ob,istatus)
            if(istatus .eq. 1)then
                write(6,*)' i4_window_ob = ',i4_window_ob
            else
                write(6,*)' Error getting i4_window_ob'
                return
            endif

        else ! append
            if(iwrite_output .ge. 1)then
                call open_lapsprd_file_append(lun_pig,i4time,ext
     1                                       ,istatus)      
                if(istatus .ne. 1)go to 888
            endif

            call get_windob_time_window('CDW',i4_window_ob,istatus)
            if(istatus .eq. 1)then
                write(6,*)' i4_window_ob = ',i4_window_ob
            else
                write(6,*)' Error getting i4_window_ob'
                return
            endif

        endif

        write(6,*)
        write(6,*)'             Reading Point Obs: ',ext_in
        write(6,*)
     1  '   n   i  j  k    u      v       dd     ff      azi    ran '

10      i_qc = 1

        if(n_point_obs .le. 500 .OR. 
     1     n_point_obs - (n_point_obs/10)*10 .eq. 9)then       
            iwrite2 = 1
        else
            iwrite2 = 0
        endif

        if(ext_in .eq. 'pin')then
            call read_acars_ob(lun_in,'wind',xlat,xlon,elev,dd,ff
     1                        ,asc9_tim_point,iwrite2,l_geoalt,l_eof)       
            if(elev .eq. 0.)i_qc = 0
        else
            call read_laps_cdw_wind(lun_in,xlat,xlon,pres_pa,dd,ff
     1                                           ,asc9_tim_point,l_eof)
        endif

        if(l_eof)goto900

        call cv_asc_i4time(asc9_tim_point,i4time_ob)

        if(abs(i4time_ob - i4time) .le. i4_window_ob)then

!           Climo/QC check
            if(abs(dd) .lt. 500. .and. abs(ff) .lt. 500. 
     1                           .and. i_qc .eq. 1)then

                call latlon_to_rlapsgrid(xlat,xlon,lat,lon,ni,nj
     1                                  ,ri,rj,istatus)
                i_grid = nint(ri)
                j_grid = nint(rj)

                if(i_grid .ge.  1 .and. j_grid .ge. 1 .and.
     1             i_grid .le. ni .and. j_grid .le. nj)then

!                   Point ob is in horizontal domain

                    if(ext_in .eq. 'pin')then
                       if(l_geoalt)then ! elev is geometric height MSL (WISDOM)
                          rk = height_to_zcoord2(elev,heights_3d
     1                        ,ni,nj,nk,i_grid,j_grid,istatus_rk)
                          if(istatus_rk .ne. 1)then
                             write(6,*)' WARNING: rejecting ACARS ',       
     1                          'apparently above top of domain ',elev
                          endif

                       else            ! ACARS elev is pressure altitude
                          if(abs(elev) .lt. 90000.)then
                             pres_mb = ztopsa(elev)
                             pres_pa = pres_mb * 100.
                             rk = rlevel_of_field(pres_pa,pres_3d
     1                                           ,ni,nj,nk
     1                                           ,i_grid,j_grid
     1                                           ,istatus_rk) 

                             if(istatus_rk .ne. 1)then
                                 write(6,*)
     1                           ' Bad status from rlevel_of_field'       
                             endif

                          else ! way too high
                             write(6,*)' WARNING: rejecting ACARS ',       
     1                          'apparently above top of domain ',elev
                             istatus_rk = 0

                          endif ! ht < 90000m

                       endif

                       weight_ob = weight_pirep

                    else ! cdw
!                      rk = zcoord_of_pressure(pres_pa)
!                      istatus_rk = 1

                       rk = rlevel_of_field(pres_pa,pres_3d
     1                                     ,ni,nj,nk
     1                                     ,i_grid,j_grid
     1                                     ,istatus_rk) 

                       if(istatus_rk .ne. 1)then
                           write(6,*)
     1                          ' Bad status from rlevel_of_field'       
                       endif

                       weight_ob = weight_cdw

                    endif

                    if (istatus_rk .eq. 1) k_grid = nint(rk)

                    if(      istatus_rk .eq. 1
     1                 .and. k_grid     .le. nk
     1                 .and. k_grid     .ge. 1 )then ! Ob is in vertical domain

                        n_point_obs = n_point_obs + 1

                        if(n_point_obs .gt. N_POINT)then
                           write(6,*)' Warning: Too many point obs, '       
     1                              ,'limit is ',N_POINT
                           istatus = 0
                           return
                        endif

                        if(n_point_obs .le. 100 .OR. 
     1                     n_point_obs .eq. (n_point_obs/100) * 100)then
                            iwrite = 1
                        else
                            iwrite = 0
                        endif

                        point_i(n_point_obs) = i_grid
                        point_j(n_point_obs) = j_grid

                        call disp_to_uv(dd,ff,u_temp,v_temp)

                        point_k(n_point_obs) = k_grid

                        call get_time_term(u_mdl_bkg_4d,ni,nj,nk
     1                                    ,NTMIN,NTMAX
     1                                    ,i_grid,j_grid,k_grid
     1                                    ,i4time,i4time_ob
     1                                    ,u_time_interp,u_diff,istatus)       

                        call get_time_term(v_mdl_bkg_4d,ni,nj,nk
     1                                    ,NTMIN,NTMAX
     1                                    ,i_grid,j_grid,k_grid
     1                                    ,i4time,i4time_ob
     1                                    ,v_time_interp,v_diff,istatus)       

!                       u_diff = du/dt * [t(ob) - t(anal)]
                        point_u(n_point_obs) = u_temp - u_diff
                        point_v(n_point_obs) = v_temp - v_diff

                        if(iwrite_output .ge. 1)then
                            if(ext_in .eq. 'pin' .and. l_geoalt)then
                                write(lun_pig,91)ri,rj,rk,dd,ff,'wis'
 91                             format(1x,5f8.1,1x,a3)
                                n_wisdom = n_wisdom + 1
                            else
                                write(lun_pig,91)ri,rj,rk,dd,ff,ext_in
                            endif
                        endif

                        if(iwrite .eq. 1)write(6,101)xlat,xlon,dd,ff,rk
     1          ,u_temp,v_temp,point_u(n_point_obs),point_v(n_point_obs)
101                     format(2f8.2,2f8.1,f8.1,4f8.2)

!                 ***   Remap point observation to LAPS observation grid

                        grid_laps_u
     1  (point_i(n_point_obs),point_j(n_point_obs),point_k(n_point_obs))      
     1  = point_u(n_point_obs)

                        grid_laps_v
     1  (point_i(n_point_obs),point_j(n_point_obs),point_k(n_point_obs))
     1  = point_v(n_point_obs)

                        grid_laps_wt
     1  (point_i(n_point_obs),point_j(n_point_obs),point_k(n_point_obs))       
     1  = weight_ob

!                       Add to data structure
                        nobs_point = nobs_point + 1
                        obs_point(nobs_point)%i = i_grid
                        obs_point(nobs_point)%j = j_grid
                        obs_point(nobs_point)%k = k_grid
                        obs_point(nobs_point)%ri = ri	! Yuanfu
                        obs_point(nobs_point)%rj = rj	! Yuanfu
                        obs_point(nobs_point)%rk = rk ! k_grid ! Yuanfu
                        obs_point(nobs_point)%valuef(1) = u_temp-u_diff       
                        obs_point(nobs_point)%valuef(2) = v_temp-v_diff       
                        obs_point(nobs_point)%weight = weight_ob
                        obs_point(nobs_point)%vert_rad_rat = 1.0
                        obs_point(nobs_point)%type   = ext_in
                        obs_point(nobs_point)%file   = ext_in
                        obs_point(nobs_point)%i4time = i4time_ob

                    else  ! Out of vertical bounds
                      iwrite = 0
                       write(6,*)' WARNING: rejecting point ob ',       
     1                        'apparently outside vertical domain ',elev  

                    endif ! In vertical bounds

                    if(iwrite .eq. 1)write(6,20)n_point_obs,
     1                 point_i(n_point_obs),
     1                 point_j(n_point_obs),
     1                 point_k(n_point_obs),
     1                 point_u(n_point_obs),
     1                 point_v(n_point_obs),
     1                 dd,ff
20                  format(i5,1x,3i4,2f7.1,2x,2f7.1,2x,2f7.1,2x,2f7.1)

                else
                    if(iwrite .eq. 1)write(6,*)
     1              ' Out of horizontal bounds',i_grid,j_grid        

                endif ! In horizontal bounds

            endif ! Good data (passed gross climo/qc check)

        else
            write(6,*)' Out of temporal bounds'
     1                              ,abs(i4time_ob - i4time)

        endif ! In temporal bounds

100     goto10

900     write(6,*)' End of RDPOINT ',ext_in,' file: '
     1           ,'Cumulative # obs = ',n_point_obs

        if(n_wisdom .gt. 0)then
            write(6,*)'# of WISDOM obs = ',n_wisdom
        endif

        close(lun_in)
        close(lun_pig)

        istatus = 1

        return

999     write(6,*)' No point ob data present for ',ext_in
        write(6,*)' End of RDPOINT ',ext_in,' file: '
     1           ,'Cumulative # obs = ',n_point_obs
        istatus = 1
        return


888     write(6,*)' Open error for PIG file'
        istatus = 0
        close(lun_pig)
        return

        end


        subroutine read_laps_cdw_wind(lun,xlat,xlon,pres,dd,ff
     1                                          ,asc9_tim_point,l_eof)

        real pres ! pa
        real dd   ! degrees (99999. is missing)
        real ff   ! meters/sec (99999. is missing)

        character*9 asc9_tim_point

        logical l_eof

        l_eof = .false.

100     read(lun,895,err=100,end=900)xlat,xlon,pres,dd,ff,asc9_tim_point       
895     FORMAT(f8.3,f10.3,f8.0,f6.0,f6.1,2x,a9)

        return

 900    l_eof = .true.

        return
        end



       subroutine get_time_term(field_4d,NX_L,NY_L,NZ_L,NTMIN,NTMAX
     1                     ,i,j,k,i4time_sys,i4time_interp
     1                     ,field_interp,field_diff
     1                     ,istatus)

!      Steve Albers 1999

!      This routine does a time interpolation

       real field_4d(NX_L,NY_L,NZ_L,NTMIN:NTMAX)

       call get_laps_cycle_time(laps_cycle_time,istatus)
       if(istatus .ne. 1)then
           write(6,*)' ERROR in get_time_term (laps_cycle_time)'
           return
       endif

       rcycles = float(i4time_interp - i4time_sys)
     1         / float(laps_cycle_time)

       if(rcycles .gt. float(NTMAX))then
           write(6,*)' Warning in get_time_term, rcycles = ',rcycles
           rcycles = NTMAX
           istatus = -1
       endif

       if(rcycles .lt. float(NTMIN))then
           write(6,*)' Warning in get_time_term, rcycles = ',rcycles
           rcycles = NTMIN
           istatus = -1
       endif

       itlow = nint(rcycles - 0.5) ! round down to nearest integer
       itlow = min(itlow,NTMAX-1)
       itlow = max(itlow,NTMIN)

       ithigh = itlow + 1

       frac_t = rcycles - float(itlow)

       field_interp = field_4d(i,j,k,itlow)  * (1.0 - frac_t)
     1              + field_4d(i,j,k,ithigh) * frac_t

!      du/dt * [t(ob) - t(anal)]
       field_diff = field_interp - field_4d(i,j,k,0) 

       istatus = 1
       return

       end
