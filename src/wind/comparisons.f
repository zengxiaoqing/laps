


        subroutine comparisons (
     1                upass1,vpass1,istat_radar_vel,max_radars,
     1                grid_ra_vel,rlat_radar,rlon_radar,rheight_radar,
     1                lat,lon,
     1                uanl,vanl,u_mdl_curr,v_mdl_curr,
     1                istat_bal,
     1                ni,nj,nk,r_missing_data,
     1                weight_pirep,weight_prof,weight_sfc,
     1                grid_laps_u,grid_laps_v,grid_laps_wt,
!    1                meso_i,meso_j,meso_k,meso_u,meso_v,N_MESO,
!    1                n_meso_obs,
!    1                sao_i,sao_j,sao_k,sao_u,sao_v,N_SAO,n_sao_obs,
     1                n_radars)

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
C          uanl
C          vanl
C          u_mdl_curr
C          v_mdl_curr
C          istat_bal
C!          ubal
C!          vbal
C          n_radars
C
C  outputs: None
C
C*********************************************************************

         implicit none

C************* Include Section ***************************************

C***************** Declarations **************************************

        integer istat_radar_vel,istat_bal
        integer l,n_radars,ni,nj,nk,max_radars

        real*4 rlat_radar(max_radars),rlon_radar(max_radars)
     1                     ,rheight_radar(max_radars)

        real*4 lat(ni,nj),lon(ni,nj)

        real*4 upass1(ni,nj,nk),vpass1(ni,nj,nk)
        real*4 grid_ra_vel(ni,nj,nk,max_radars),r_missing_data
        real*4 weight_pirep,weight_prof,weight_sfc,weight_cdw
        real*4 grid_laps_u(ni,nj,nk),grid_laps_v(ni,nj,nk)
     1                                          ,grid_laps_wt(ni,nj,nk)

        real*4 uanl(ni,nj,nk),vanl(ni,nj,nk)
        real*4 u_mdl_curr(ni,nj,nk),v_mdl_curr(ni,nj,nk)

        real*4 rms_fg_cdw,rms_fg_sfc,rms_laps_cdw,rms_laps_sfc
        real*4 rms_fg_vr,rms_laps_vr,rms_laps_prof,rms_fg_pirep
        real*4 rms_laps_pirep,rms_maps_prof,rms_maps_pirep,rms_fg_laps
        real*4 rms_laps_maps,rms_maps_sfc

C********************************************************************

        write(6,*)
        write(6,*)'  Comparing LAPS Analysis to CDW Obs (passing QC)'       
        call comp_grid_windobs(uanl,vanl,ni,nj,nk
     1        ,grid_laps_u,grid_laps_v,grid_laps_wt,weight_cdw
     1        ,' FG ','CDW ',r_missing_data,rms_fg_cdw)

        write(6,*)
        write(6,*)'  Comparing LAPS First Pass to SFC Obs (passing QC)'       
        call comp_laps_sfc(upass1,vpass1,ni,nj,nk
     1        ,grid_laps_u,grid_laps_v,grid_laps_wt,weight_sfc
     1        ,r_missing_data,rms_fg_sfc)

        write(6,*)
        write(6,*)'  Comparing LAPS Analysis to CDW Obs (passing QC)'       
        call comp_grid_windobs(uanl,vanl,ni,nj,nk
     1        ,grid_laps_u,grid_laps_v,grid_laps_wt,weight_cdw
     1        ,'LAPS','CDW ',r_missing_data,rms_laps_cdw)

        write(6,*)
        write(6,*)'  Comparing LAPS Analysis to SFC Obs (passing QC)'       
        call comp_laps_sfc(uanl,vanl,ni,nj,nk
     1        ,grid_laps_u,grid_laps_v,grid_laps_wt,weight_sfc
     1        ,r_missing_data,rms_laps_sfc)

        do l = 1,n_radars

            write(6,*)
            write(6,*)'  Comparing LAPS First Pass to Radial Velocities'       
     1                                                    ,' Radar #',l
            if(istat_radar_vel .eq. 1)
     1        call comp_laps_vr(grid_ra_vel(1,1,1,l),upass1,vpass1
     1          ,ni,nj,nk,r_missing_data,rms_fg_vr
     1          ,lat,lon,rlat_radar(l),rlon_radar(l),rheight_radar(l))

            write(6,*)
            write(6,*)'  Comparing LAPS Analysis to Radial Velocities'
     1                                                    ,' Radar #',l
            if(istat_radar_vel .eq. 1)
     1        call comp_laps_vr(grid_ra_vel(1,1,1,l),uanl,vanl,ni,nj,nk       
     1          ,r_missing_data,rms_laps_vr
     1          ,lat,lon,rlat_radar(l),rlon_radar(l),rheight_radar(l))


        enddo ! l

        write(6,*)
        write(6,*)'  Comparing LAPS First Pass to Profiler'
        call comp_laps_prof(upass1,vpass1,ni,nj,nk
     1        ,r_missing_data,weight_prof,grid_laps_u,grid_laps_v
     1        ,grid_laps_wt,rms_laps_prof)

        write(6,*)
        write(6,*)'  Comparing LAPS Analysis to Profiler'
        call comp_laps_prof(uanl,vanl,ni,nj,nk
     1        ,r_missing_data,weight_prof,grid_laps_u,grid_laps_v
     1        ,grid_laps_wt,rms_laps_prof)

        write(6,*)
        write(6,*)'  Comparing LAPS First Pass to Pireps'
        call comp_laps_pirep(upass1,vpass1,ni,nj,nk
     1        ,r_missing_data,weight_pirep,grid_laps_u,grid_laps_v
     1        ,grid_laps_wt,rms_fg_pirep)

        write(6,*)
        write(6,*)'  Comparing LAPS Analysis to Pireps'
        call comp_laps_pirep(uanl,vanl,ni,nj,nk,r_missing_data
     1         ,weight_pirep,grid_laps_u,grid_laps_v,grid_laps_wt
     1         ,rms_laps_pirep)

        write(6,*)
        write(6,*)'  Comparing Model Background to Profiler'
        call comp_maps_prof(u_mdl_curr,v_mdl_curr,ni,nj,nk
     1        ,grid_laps_u,grid_laps_v,grid_laps_wt,weight_prof
     1        ,r_missing_data,rms_maps_prof)

        write(6,*)
        write(6,*)'  Comparing Model Background to SFC Obs'
        call comp_maps_sao(u_mdl_curr,v_mdl_curr,ni,nj,nk
     1        ,grid_laps_u,grid_laps_v,grid_laps_wt,weight_sfc
     1        ,r_missing_data,rms_maps_sfc)

        write(6,*)
        write(6,*)'  Comparing Model Background Analysis to Pireps'
        call comp_maps_pirep(u_mdl_curr,v_mdl_curr,ni,nj,nk
     1        ,grid_laps_u,grid_laps_v,grid_laps_wt,weight_pirep
     1        ,r_missing_data,rms_maps_pirep)

        write(6,*)
        write(6,*)'  Comparing LAPS First Pass & Analysis'
        call comp_laps1_laps2(upass1,vpass1,uanl,vanl
     1                                  ,ni,nj,nk,rms_fg_laps)

        write(6,*)
        write(6,*)'  Comparing LAPS Analysis & MODEL'
        call comp_laps_maps(uanl,vanl,u_mdl_curr,v_mdl_curr,ni,nj,nk
     1             ,r_missing_data,rms_laps_maps)


        return
        end



        subroutine comp_grid_windobs(u_3d,v_3d,ni,nj,nk
     1  ,grid_laps_u,grid_laps_v,grid_laps_wt,weight_ob
     1  ,c_grid,c_ob,r_missing_data,rms)

        real*4 grid_laps_u(ni,nj,nk),grid_laps_v(ni,nj,nk)
        real*4 grid_laps_wt(ni,nj,nk)
        dimension u_3d(ni,nj,nk),v_3d(ni,nj,nk) 

        character*(*)c_grid,c_ob

        nobs = 0
        residualu = 0.
        residualv = 0.

        write(6,*)'Comparing ',c_ob,' Wind Obs (passing QC) to '
     1                       ,c_grid,' Grid'
        write(6,2)c_ob,c_grid
2       format(1x,'   i   j   k      ',a,'                 '
     1                                ,a,' Analysis diff')

        do jl = 1,nj

          do il = 1,ni

              do k = 1,nk

                if(grid_laps_wt(il,jl,k) .eq. weight_ob
     1          .and. u_3d(il,jl,k) .ne. r_missing_data )then
                  nobs = nobs + 1

                  diffu = u_3d(il,jl,k) - grid_laps_u(il,jl,k)
                  diffv = v_3d(il,jl,k) - grid_laps_v(il,jl,k)
                  residualu = residualu + diffu ** 2
                  residualv = residualv + diffv ** 2
                  write(6,101)il,jl,k
     1              ,grid_laps_u(il,jl,k),grid_laps_v(il,jl,k)
     1              ,u_3d(il,jl,k),v_3d(il,jl,k)
     1              ,-diffu,-diffv
101               format(1x,3i4,3(2x,2f7.1))

                endif

              enddo ! k

          enddo ! j
        enddo ! i


        if(nobs .gt. 0)then
            rmsu = sqrt(residualu/nobs)
            rmsv = sqrt(residualv/nobs)
        else
            rmsu = 0.
            rmsv = 0.
        endif

        rms  = sqrt(rmsu**2 + rmsv**2)

        write(6,102)c_ob,c_grid,nobs,rmsu,rmsv,rms
102     format(' RMS between ',a,' & ',a,' (n,rmsu,rmsv,rms) = ',
     1     i4,3f5.1)

        return

        end

