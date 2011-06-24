          
        subroutine verif_fcst_pt_3d(i4time_sys,a9time,laps_cycle_time,
     1                  ni,nj,
     1                  nk,
     1                  maxsta,
     1                  r_missing_data,
     1                  model_cycle_time_sec,
     1                  n_fcst_times,
     1                  j_status)

        include 'barnesob.inc'
        include 'windparms.inc'

        real var_anal_3d(ni,nj,nk)
        real var_fcst_3d(ni,nj,nk)
        real lat(ni,nj)
        real lon(ni,nj)
        real topo(ni,nj)
        real rlaps_land_frac(ni,nj)

        parameter (NTMIN=-1) ! read from namelist file?
        parameter (NTMAX=+1)

!       For call to get_wind_3d_obs
        real u_mdl_bkg_4d(ni,nj,nk,NTMIN:NTMAX)
        real v_mdl_bkg_4d(ni,nj,nk,NTMIN:NTMAX)
        real grid_laps_u(ni,nj,nk)
        real grid_laps_v(ni,nj,nk)
        real grid_laps_wt(ni,nj,nk)
        real heights_3d(ni,nj,nk)
        real heights_1d(nk)
        type (barnesob) :: obs_point(maxsta)   

        real k_to_f

        integer       maxbgmodels
        parameter     (maxbgmodels=10)
        character*30  c_fdda_mdl_src(maxbgmodels)

        integer max_fcst_times
        parameter (max_fcst_times=200)

        integer max_regions
        parameter (max_regions=10)

        character EXT*31, directory*255, c_model*30

        character*10  units_2d
        character*125 comment_2d
        character*3 var_2d
        character*9 a9time,a9time_valid,a9time_init
        character*24 atime_s
        character*150 hist_dir, verif_dir
        character*150 hist_file

        integer n_fields
        parameter (n_fields=2)
        character*10 ext_anal_a(n_fields), ext_fcst_a(n_fields)
        character*10 var_a(n_fields)
        character*2 c2_region

!       Specify what is being verified
        data ext_fcst_a /'fua','fua'/ ! 3-D
!       data ext_anal_a /'lps'/ ! 3-D reflectivity
        data var_a      /'U3','V3'/ ! 3-D Wind

        integer contable(0:1,0:1)

        integer maxthr
        parameter (maxthr=5)

!       Declarations from compare_analysis_to_rad
        real                   
     1                        cvr_max(ni,nj),cvr_sao_max(ni,nj)
     1        ,dbz_max_2d(ni,nj),solar_alt(ni,nj)                   

!       How much the solar radiation varies with changes in cloud fraction
        real cvr_scl_a(ni,nj) 

        real rad_clr(ni,nj)

        real airmass_2d(ni,nj)
        real pw_2d(ni,nj)
        real trans_h2o_2d(ni,nj)

        real dum_2d(ni,nj)

        real cld_snd(max_cld_snd,KCLOUD)
        integer ista_snd(max_cld_snd)
        real cld_hts(KCLOUD)

        character c_stations(maxsta)*3
        character stn(maxsta)*20

        real lat_s(maxsta), lon_s(maxsta), elev_s(maxsta)
        real var_s(maxsta)
        real rad2_s(maxsta)
        real var_fcst_s(maxsta)
        real cvr_s(maxsta)

        real dum_s(maxsta)

        real ea(maxsta)

        integer ii_s(maxsta)
        integer jj_s(maxsta)
        integer doy

        character*3 c3_discrep
        character*1 c1_c
        character title*60, ver_file*256, filename*9

!       End Declarations

        write(6,*)' Start subroutine verif_fcst_pt_3d...'
        write(6,*)' i4time_sys a9time = ',i4time_sys,a9time

        call get_laps_domain_95(ni,nj,lat,lon,topo
     1           ,rlaps_land_frac,grid_spacing_cen_m,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error getting LAPS domain'
            return
        endif

        u_mdl_bkg_4d = 0 ! Sets time term to zero when we call get_wind_3d_obs
        v_mdl_bkg_4d = 0

!       n_fcst_times = 2 ! 38

        thresh_var = 20. ! lowest threshold for this variable

        i4_initial = i4time_sys

        lun_in = 21

        rmiss = -999.

!       Get fdda_model_source from static file
        call get_fdda_model_source(c_fdda_mdl_src,n_fdda_models,istatus)

        write(6,*)' n_fdda_models = ',n_fdda_models
        write(6,*)' c_fdda_mdl_src = '
     1            ,(c_fdda_mdl_src(m),m=1,n_fdda_models)

!       if(n_fdda_models .ne. n_models + 1)then
!           write(6,*)' ERROR n_models differs from n_fdda_models '
!    1                       ,n_models,n_fdda_models
!           stop
!       endif

        do ifield = 1,n_fields

          var_2d = var_a(ifield)
          call s_len(var_2d,lenvar)

          write(6,*) 
          write(6,*) 
          write(6,*)' Processing field ',trim(var_2d) 

          call get_directory('verif',verif_dir,len_verif)

          istart = 0

          lun_out = 39

          hist_dir = verif_dir(1:len_verif)//var_2d(1:lenvar)
     1                                       //'/pt'
!    1                                       //c_model(1:len_model)
          len_hist = len_verif + 3 + lenvar ! + len_model

          call make_fnam_lp(i4_initial,a9time_init,istatus)

          hist_file = hist_dir(1:len_hist)//'/'//a9time_init
     1                                    //'.stats'     

!         if(ifield .gt. 1)then
!           hist_file = trim(hist_file)//'_'//var_a(ifield)
!         endif

          call s_len(hist_file,len_histf)
          write(6,*)'hist_file = ',hist_file(1:len_histf)
          open(lun_out,file=hist_file(1:len_histf),status='unknown')

          do imodel = 1,n_fdda_models ! move this inside the time loop

            c_model = c_fdda_mdl_src(imodel)

            if(c_model(1:3) .ne. 'lga')then

              write(6,*)
              write(6,*)' Processing model ',c_model

              call s_len(c_model,len_model)

              if(istart .eq. 1)then
                  call cv_i4tim_asc_lp(i4_initial,atime_s,istatus)
                  write(lun_out,710)atime_s
              endif

              do itime_fcst = istart,n_fcst_times

                itime = itime_fcst + 1

                i4_valid = i4_initial 
     1                   + itime_fcst * model_cycle_time_sec

                call make_fnam_lp(i4_valid,a9time_valid,istatus)

                write(6,*) 
                write(6,*)' Processing time ',a9time_valid

                write(c2_region,1)iregion
 1              format(i2.2)

                N_SAO = 100000      ! read from wind.nl or nest7grid.parms?
                N_PIREP = 100000    ! read from wind.nl or nest7grid.parms?
                MAX_PR = 3000       ! read from wind.nl?
                MAX_PR_LEVELS = 300 ! read from wind.nl?

                nobs_point = 0      ! initialize
                var_s = r_missing_data
                var_fcst_s = r_missing_data
                max_wind_obs = maxsta

!               Read wind obs        
                if(var_2d .eq. 'U3' .or. var_2d .eq. 'V3')then
                  write(6,*)' Reading analyzed height field'
                  call get_laps_3d(i4_valid
     1                          ,ni,nj,nk
     1                          ,'lt1'    ,'HT '  
     1                          ,units_2d,comment_2d
     1                          ,heights_3d 
     1                          ,istatus)
                  if(istatus .ne. 1)then
                       write(6,*)' Error reading analyzed height '
                       goto 990
                  endif
                  
                  write(6,*)' Reading 3D wind obs - i4_valid = '
     1                                             ,i4_valid
                  call get_wind_3d_obs(
     1            ni,nj,nk,                                       ! I
     1            r_missing_data,i2_missing_data,                 ! I
     1            i4_valid,heights_3d,heights_1d,                 ! I
     1            MAX_PR,MAX_PR_LEVELS,weight_prof,l_use_raob,    ! I
     1            l_use_cdw,                                      ! I
     1            N_SAO,N_PIREP,                                  ! I
     1            lat,lon,                                        ! I
     1            NTMIN,NTMAX,                                    ! I
     1            u_mdl_bkg_4d, v_mdl_bkg_4d,                     ! I
     1            grid_laps_u,grid_laps_v,grid_laps_wt,           ! O
     1            max_wind_obs,obs_point,nobs_point,              ! I/O
     1            rlat_radar,rlon_radar,rheight_radar,            ! I
     1            istat_radar_vel,n_vel_grids,                    ! I
     1            istatus_remap_pro,                              ! O
     1            istatus                )                        ! O

                endif

                write(6,*)' number of obs ',nobs_point      

                if(var_2d .eq. 'SWI')then
                  var_s = solar_s
                  threshval = 0.
                elseif(var_2d .eq. 'TSF')then
                  var_s = t_s
                  threshval = -99.9
                endif

                if(.true.)then

                  write(6,*)' Reading forecast field'

!                 Read forecast field
                  ext = ext_fcst_a(ifield)
                  call get_directory(ext,directory,len_dir)

                  DIRECTORY=directory(1:len_dir)//c_model(1:len_model)
     1                                          //'/'

                  call get_lapsdata_3d(i4_initial,i4_valid
     1                          ,ni,nj,nk
     1                          ,DIRECTORY,var_2d
     1                          ,units_2d,comment_2d
     1                          ,var_fcst_3d
     1                          ,istatus)
                  if(istatus .ne. 1)then
                       write(6,*)' Error reading 3D Forecast for '
     1                           ,var_2d
                       xbar = rmiss
                       ybar = rmiss
                       std = rmiss
                       goto 900
                  endif

                endif ! .true.

                cnt = 0.

                do ista = 1,nobs_point
                  if(var_2d .eq. 'V3')then
                    var_s(ista) = obs_point(ista)%valuef(2)
                  else
                    var_s(ista) = obs_point(ista)%valuef(1)
                  endif

                  var_fcst_s(ista) = r_missing_data 
                  rad2_s(ista) = r_missing_data 

                  write(6,*)ista,var_s(ista)

                  if(var_s(ista) .ne. r_missing_data)then
                    call latlon_to_rlapsgrid(lat_s(ista),lon_s(ista),lat
     1                          ,lon,ni,nj,ri,rj,istatus)

                    i_i = obs_point(ista)%i
                    i_j = obs_point(ista)%j
                    i_k = obs_point(ista)%k

                    ii_s(ista) = i_i
                    jj_s(ista) = i_j

                    if(i_i .ge. 3 .and. i_i .le. ni-2 .and.
     1                 i_j .ge. 3 .and. i_j .le. nj-2            )then

                      if(iwrite .eq. iwrite/20*20)then
                        write(6,*)'sv '
                        write(6,*)'sv Sta   i    j   VIS frac tb8_k  '
     1                  //'t_gnd_k t_sfc_k cv_s_mx cvr_mx '
     1                  //'solalt 9pt  rad_fc '
     1                  //'rad_ob rad_th ratio cv_sol  df'
                      endif

                      iwrite = iwrite + 1

                      if(var_2d .eq. 'U3' .or. var_2d .eq. 'V3')then
                        var_fcst_s(ista) = var_fcst_3d(i_i,i_j,i_k)
                      endif
 
                      sumobs = sumobs + var_s(ista)
                      sumanl = sumanl + var_fcst_s(ista)
                      cnt = cnt + 1.

                      cvr_s(ista) = cvr_max(i_i,i_j)

1112                endif ! ob is in domain
                  endif ! ista .ne. 0 (valid value)
                enddo ! ista

                write(6,*)
                write(6,*)' Generic stats, cnt = ',nint(cnt)
                if(var_2d .eq. 'U3')then
                  call stats_1d(maxsta,var_fcst_s,var_s 
     1                   ,'U Wind Component:       '
     1                   ,a_t,b_t,xbar,ybar
     1                   ,bias,std,r_missing_data,istatus)
                elseif(var_2d .eq. 'V3')then
                  call stats_1d(maxsta,var_fcst_s,var_s
     1                   ,'V Wind Component:       '
     1                   ,a_t,b_t,xbar,ybar
     1                   ,bias,std,r_missing_data,istatus)
                endif

900             write(6,*)
                write(6,*)' Writing to lun_out ',lun_out
                call cv_i4tim_asc_lp(i4_valid,atime_s,istatus)
                write(6,710)atime_s,xbar,ybar,std
                write(lun_out,710)atime_s,xbar,ybar,std
!               write(39,710)atime_s,xbar,ybar,std
!               write(38,710)atime_s,xbar,ybar,std
710             format(1x,a24,3f10.3)

990             continue

              enddo ! itime_fcst

              write(lun_out,*)
              write(lun_out,*)

            endif ! c_model .ne. lga

          enddo ! model

          write(6,*)' Closing lun_out ',lun_out
          close (lun_out) 

        enddo ! fields

 999    write(6,*)' End of subroutine verif_fcst_pt_3d'

        return

        end
