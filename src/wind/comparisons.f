comparisons.f.cdis   
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



        subroutine compare_wind (
     1                upass1,vpass1,cgrid,
     1                istat_radar_vel,max_radars,grid_ra_vel,n_radars,
     1                rlat_radar,rlon_radar,rheight_radar,
     1                lat,lon,
     1                ni,nj,nk,r_missing_data,
     1                obs_barnes,max_obs,ncnt_total,l_point_struct,
     1                weight_pirep,weight_prof,weight_sfc,weight_cdw,       
     1                grid_laps_u,grid_laps_v,grid_laps_wt,istatus)

C****************************************************************************
C
C  Purpose: Provide a single point out of lapswind_anal to call
C           diagnostic comparision routines.
C
C
C  Inputs: upass1
C          vpass1
C          istat_radar_vel
C          grid_ra_vel
C          rlat_radar
C          rlon_radar
C          rheight_radar
C          n_radars
C
C  outputs: None
C
C*********************************************************************

C***************** Declarations **************************************
        include 'barnesob.inc'
        type (barnesob) obs_barnes(max_obs)      

        integer istat_radar_vel
        integer l,n_radars,ni,nj,nk,max_radars

        real*4 rlat_radar(max_radars),rlon_radar(max_radars)
     1                     ,rheight_radar(max_radars)

        real*4 lat(ni,nj),lon(ni,nj)

        real*4 upass1(ni,nj,nk),vpass1(ni,nj,nk)
        real*4 grid_ra_vel(ni,nj,nk,max_radars),r_missing_data
        real*4 weight_pirep,weight_prof,weight_sfc,weight_cdw
        real*4 grid_laps_u(ni,nj,nk),grid_laps_v(ni,nj,nk)
     1                                          ,grid_laps_wt(ni,nj,nk)

        character*4 cgrid
        logical l_parse, l_point_struct

C********************************************************************

        write(6,*)' Subroutine compare_wind...',cgrid

        write(6,*)
        if(l_parse(cgrid,'FG'))then
            write(6,*)'  Comparing ',cgrid,' to SFC Obs (prior to QC)'       
        else
            write(6,*)'  Comparing ',cgrid,' to SFC Obs (passing QC)'       
        endif
        call comp_grid_windobs(upass1,vpass1,ni,nj,nk
     1        ,grid_laps_u,grid_laps_v,grid_laps_wt,weight_sfc
     1        ,obs_barnes,max_obs,ncnt_total,l_point_struct
     1        ,cgrid,'SFC ',r_missing_data,rms)

        write(6,*)
        write(6,*)'  Comparing ',cgrid,' to Profiler'       
        call comp_grid_windobs(upass1,vpass1,ni,nj,nk
     1        ,grid_laps_u,grid_laps_v,grid_laps_wt,weight_prof
     1        ,obs_barnes,max_obs,ncnt_total,l_point_struct
     1        ,cgrid,'PROF',r_missing_data,rms)

        write(6,*)
        write(6,*)'  Comparing ',cgrid,' to Pireps'       
        call comp_grid_windobs(upass1,vpass1,ni,nj,nk
     1        ,grid_laps_u,grid_laps_v,grid_laps_wt,weight_pirep
     1        ,obs_barnes,max_obs,ncnt_total,l_point_struct
     1        ,cgrid,'PIN ',r_missing_data,rms)

        write(6,*)
        if(l_parse(cgrid,'FG'))then
            write(6,*)'  Comparing ',cgrid,' to CDW Obs (prior to QC)'       
        else
            write(6,*)'  Comparing ',cgrid,' to CDW Obs (passing QC)'       
        endif
        call comp_grid_windobs(upass1,vpass1,ni,nj,nk
     1        ,grid_laps_u,grid_laps_v,grid_laps_wt,weight_cdw
     1        ,obs_barnes,max_obs,ncnt_total,l_point_struct
     1        ,cgrid,'CDW ',r_missing_data,rms)

!       write(6,*)
!       write(6,*)'  Comparing LAPS First Pass & Analysis'
!       call comp_laps1_laps2(upass1,vpass1,uanl,vanl
!    1                                  ,ni,nj,nk,rms_fg_laps)

!       write(6,*)
!       write(6,*)'  Comparing LAPS Analysis & MODEL'
!       call comp_laps_maps(uanl,vanl,u_mdl_curr,v_mdl_curr,ni,nj,nk
!    1             ,r_missing_data,rms_laps_maps)

        do l = 1,n_radars
            write(6,*)
            write(6,*)'  Comparing ',cgrid,' to Radial Velocities'       
     1                                                    ,' Radar #',l
            if(istat_radar_vel .eq. 1)
     1        call comp_laps_vr(grid_ra_vel(1,1,1,l),upass1,vpass1
     1          ,ni,nj,nk,r_missing_data,cgrid,rms_fg_vr
     1          ,lat,lon,rlat_radar(l),rlon_radar(l),rheight_radar(l))

        enddo ! l

        return
        end



        subroutine comp_grid_windobs(u_3d,v_3d,ni,nj,nk
     1  ,grid_laps_u,grid_laps_v,grid_laps_wt,weight_ob
     1  ,obs_barnes,max_obs,ncnt_total,l_point_struct
     1  ,c_grid,c_ob,r_missing_data,rms)

        include 'barnesob.inc'
        type (barnesob) obs_barnes(max_obs)      

        real*4 grid_laps_u(ni,nj,nk),grid_laps_v(ni,nj,nk)
        real*4 grid_laps_wt(ni,nj,nk)
        dimension u_3d(ni,nj,nk),v_3d(ni,nj,nk) 

        character*(*)c_grid,c_ob
        character*12 c_ob_loc

        logical l_point_struct

        nobs = 0
        residualu = 0.
        residualv = 0.
        sumu = 0.
        sumv = 0.
        sumsp = 0.

        write(6,*)'Comparing ',c_ob,' Wind Obs (passing QC) to '
     1                       ,c_grid,' Grid'
        write(6,2)c_ob,c_grid
2       format(1x,'   i   j   k      ',a,' Ob          '
     1                                ,a,' Analysis     diff')

        if(l_point_struct)then
           c_ob_loc = c_ob
           call downcase(c_ob_loc,c_ob_loc)
           call left_justify(c_ob_loc)
 
           do iob = 1,ncnt_total
              if(obs_barnes(iob)%type .eq. c_ob_loc)then
                  nobs = nobs + 1
                  il = obs_barnes(iob)%i
                  jl = obs_barnes(iob)%j
                  k = obs_barnes(iob)%k
                  diffu = obs_barnes(iob)%value(1) - u_3d(il,jl,k) 
                  diffv = obs_barnes(iob)%value(2) - v_3d(il,jl,k) 

                  sumu = sumu + diffu
                  sumv = sumv + diffv

                  residualu = residualu + diffu ** 2
                  residualv = residualv + diffv ** 2

!                 Calculate Speed Difference (Ob - Grid)
                  call uv_to_disp(obs_barnes(iob)%value(1),
     1                            obs_barnes(iob)%value(2),
     1                            di_dum,
     1                            grid_laps_sp)

                  call uv_to_disp(u_3d(il,jl,k),
     1                            v_3d(il,jl,k),
     1                            di_dum,
     1                            sp_3d)

                  sumsp = sumsp + (grid_laps_sp - sp_3d)
                  if(nobs .le. 200 .OR. nobs .eq. (nobs/10)*10)then
                      write(6,101)il,jl,k
     1                ,obs_barnes(iob)%value(1),obs_barnes(iob)%value(2)       
     1                ,u_3d(il,jl,k),v_3d(il,jl,k)
     1                ,diffu,diffv
101                   format(1x,3i4,3(2x,2f7.1))
                  endif

              endif ! obstype match

           enddo ! iob

        else
           do jl = 1,nj

           do il = 1,ni

              do k = 1,nk

                if(grid_laps_wt(il,jl,k) .eq. weight_ob
     1             .and. u_3d(il,jl,k) .ne. r_missing_data )then
                  nobs = nobs + 1

!                 Calculate U/V Difference (Ob - Grid)
                  diffu = grid_laps_u(il,jl,k) - u_3d(il,jl,k) 
                  diffv = grid_laps_v(il,jl,k) - v_3d(il,jl,k) 

                  sumu = sumu + diffu
                  sumv = sumv + diffv

                  residualu = residualu + diffu ** 2
                  residualv = residualv + diffv ** 2

!                 Calculate Speed Difference (Ob - Grid)
                  call uv_to_disp(grid_laps_u(il,jl,k),
     1                            grid_laps_v(il,jl,k),
     1                            di_dum,
     1                            grid_laps_sp)

                  call uv_to_disp(u_3d(il,jl,k),
     1                            v_3d(il,jl,k),
     1                            di_dum,
     1                            sp_3d)

                  sumsp = sumsp + (grid_laps_sp - sp_3d)

                  if(nobs .le. 200 .OR. nobs .eq. (nobs/10)*10)then
                      write(6,101)il,jl,k
     1                  ,grid_laps_u(il,jl,k),grid_laps_v(il,jl,k)
     1                  ,u_3d(il,jl,k),v_3d(il,jl,k)
     1                  ,diffu,diffv
                  endif

                endif

              enddo ! k

           enddo ! il
           enddo ! jl

        endif ! l_point_struct

        if(nobs .gt. 0)then
            rmsu = sqrt(residualu/nobs)
            rmsv = sqrt(residualv/nobs)
            biasu = sumu / float(nobs)
            biasv = sumv / float(nobs)
            biassp = sumsp / float(nobs)

        else
            rmsu = 0.
            rmsv = 0.
            biasu = 0.
            biasv = 0.
            biassp = 0.

        endif

        rms  = sqrt(rmsu**2 + rmsv**2)

        write(6,102)c_ob,c_grid,nobs,biasu,biasv,biassp,rmsu,rmsv,rms
102     format(' BIAS/RMS between '
     1        ,a,' & ',a,' (n,biasu,biasv,biassp,rmsu,rmsv,rms) = '
     1        ,i4,6f5.1)

        return

        end

