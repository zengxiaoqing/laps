          
        subroutine verif_fcst_pt_3d(i4time_sys,a9time,laps_cycle_time,
     1                  ni,nj,
     1                  nk,
     1                  maxsta,max_obs,
     1                  r_missing_data,
     1                  model_verif_intvl,
     1                  n_fcst_times,
     1                  j_status)

        include 'barnesob.inc'
        include 'windparms.inc'
        include 'tempobs.inc'   

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
        real pres_3d(ni,nj,nk)
        type (barnesob) :: obs_point(maxsta)   

        real k_to_c

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

        integer n_vfields
        parameter (n_vfields=4)
        character*10 ext_anal_a(n_vfields), ext_fcst_a(n_vfields)
        character*10 var_a(n_vfields)

!       Specify what is being verified
        data ext_fcst_a /'fua','fua','   ','fua'/ ! 3-D
!       data ext_anal_a /'lps'/ ! 3-D reflectivity
        data var_a      /'U3','V3','W3','T3'/               

        real rms_a (n_vfields,maxbgmodels,0:max_fcst_times)

        integer contable(0:1,0:1)

        integer maxthr
        parameter (maxthr=5)

!       Declarations from compare_analysis_to_rad
        real                   
     1                        cvr_max(ni,nj)                    
     1        ,dbz_max_2d(ni,nj),solar_alt(ni,nj)                   

!       How much the solar radiation varies with changes in cloud fraction
        real rad_clr(ni,nj)

        real airmass_2d(ni,nj)
        real pw_2d(ni,nj)
        real trans_h2o_2d(ni,nj)

        real dum_2d(ni,nj)

        character c_stations(maxsta)*3

        real lat_s(maxsta), lon_s(maxsta), elev_s(maxsta)
        real var_s(maxsta)
        real var_fcst_s(maxsta)

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

        iwrite = 0

        u_mdl_bkg_4d = 0 ! Sets time term to zero when we call get_wind_3d_obs
        v_mdl_bkg_4d = 0

!       n_fcst_times = 2 ! 38

        thresh_var = 20. ! lowest threshold for this variable

        i4_initial = i4time_sys

        lun_in = 21

        rmiss = -99.9

        rms_a = r_missing_data

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

!       do ifield = 1,3                  
        do ifield = 1,n_vfields

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

!         Write comment with model member names
          write(lun_out,902)(trim(c_fdda_mdl_src(imodel))
     1                     ,imodel=2,n_fdda_models)
902       format('# ',30(1x,a))

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
     1                   + itime_fcst * model_verif_intvl

                call make_fnam_lp(i4_valid,a9time_valid,istatus)

                write(6,*) 
                write(6,*)' Processing time ',a9time_valid

                if(trim(var_2d) .eq. 'W3')then
!                   call cv_i4tim_asc_lp(i4_valid,atime_s,istatus)
                    goto 1200
                endif

                N_SAO = 100000      ! read from wind.nl or nest7grid.parms?
                N_PIREP = 100000    ! read from wind.nl or nest7grid.parms?
                MAX_PR = 3000       ! read from wind.nl?
                MAX_PR_LEVELS = 300 ! read from wind.nl?

                nobs_point = 0      ! initialize
                var_s = r_missing_data
                var_fcst_s = r_missing_data
                max_wind_obs = maxsta

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

!               Read obs        
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
                       write(6,*)
     1                 ' It will try reading forecasted height '
   
                       call get_lapsdata_3d(i4_initial,i4_valid
     1                 ,ni,nj,nk,DIRECTORY,'HT ',units_2d,comment_2d
     1                 ,heights_3d,istatus)
                       if(istatus .ne. 1)then
                            write(6,*)' Error reading fcst height '
                            goto 990
                       endif
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
     1            lat,lon,topo,                                   ! I
     1            NTMIN,NTMAX,                                    ! I
     1            u_mdl_bkg_4d, v_mdl_bkg_4d,                     ! I
     1            grid_laps_u,grid_laps_v,grid_laps_wt,           ! O
     1            max_wind_obs,obs_point,nobs_point,              ! I/O
     1            rlat_radar,rlon_radar,rheight_radar,            ! I
     1            istat_radar_vel,n_vel_grids,                    ! I
     1            istatus_remap_pro,                              ! O
     1            istatus                )                        ! O

                elseif(var_2d .eq. 'T3')then ! read temperature obs
                  write(6,*)' Reading analyzed height field'
                  call get_laps_3d(i4_valid
     1                          ,ni,nj,nk
     1                          ,'lt1'    ,'HT '  
     1                          ,units_2d,comment_2d
     1                          ,heights_3d 
     1                          ,istatus)
                  if(istatus .ne. 1)then
                       write(6,*)' Error reading analyzed height '
                       write(6,*)
     1                 ' It will try reading forecasted height '

                       call get_lapsdata_3d(i4_initial,i4_valid
     1                 ,ni,nj,nk,DIRECTORY,'HT ',units_2d,comment_2d
     1                 ,heights_3d,istatus)
                       if(istatus .ne. 1)then
                            write(6,*)' Error reading fcst height '
                            goto 990
                       endif
                  endif

                  call get_pres_3d(i4_valid,ni,nj,nk,pres_3d,istatus)
                  if(istatus .ne. 1)then
                    write(6,*)
     1                ' Error: Bad status returned from get_pres_3d'
                    goto 990
                  endif

                  call get_meso_sao_pirep(MAX_SFC,dum,MAX_ACARS,istatus)                 

                  nobs_point = 0 ! initialize

                  write(6,*)' Reading 3D ACARS temp obs - i4_valid = '
     1                                                   ,i4_valid
                  call rd_acars_t(i4_valid,heights_3d,var_fcst_3d   ! I
     1                       ,pres_3d                               ! I
     1                       ,MAX_ACARS                             ! I
     1                       ,n_good_acars                          ! O
     1                       ,'pin'                                 ! I
     1                       ,ni,nj,nk                              ! I
     1                       ,lat,lon,r_missing_data                ! I
     1                       ,temp_obs,max_obs,nobs_point           ! I/O
     1                       ,istatus)                              ! O
                  if(istatus .ne. 1)then
                    write(6,*)' bad status return from rd_acars_t'
                    return
                  endif

                endif

                write(6,*)' number of obs ',nobs_point      

                cnt = 0.

                do ista = 1,nobs_point
                  if(var_2d .eq. 'U3')then
                    var_s(ista) = obs_point(ista)%valuef(1)
                  elseif(var_2d .eq. 'V3')then
                    var_s(ista) = obs_point(ista)%valuef(2)
                  elseif(var_2d .eq. 'T3')then
                    var_s(ista) = temp_obs(ista,i_ob_grid)             
                  endif

                  var_fcst_s(ista) = r_missing_data 

                  if(ista .le. 100 .or. ista .eq. (ista/100)*100)then
                      write(6,*)ista,var_s(ista)
                  endif

                  if(var_s(ista) .ne. r_missing_data)then
                    call latlon_to_rlapsgrid(lat_s(ista),lon_s(ista),lat
     1                          ,lon,ni,nj,ri,rj,istatus)

                    if(var_2d .eq. 'T3')then
                      i_g = temp_obs(ista,i_i)             
                      j_g = temp_obs(ista,i_j)
                      k_g = temp_obs(ista,i_k)
                    else
                      i_g = obs_point(ista)%i
                      j_g = obs_point(ista)%j
                      k_g = obs_point(ista)%k
                    endif

                    ii_s(ista) = i_g
                    jj_s(ista) = j_g 

                    if(i_g .ge. 3 .and. i_g .le. ni-2 .and.
     1                 j_g .ge. 3 .and. j_g .le. nj-2            )then

                      if(k_g .le. nk)then

!                       if(var_2d .eq. 'U3' .or. var_2d .eq. 'V3')then
                          var_fcst_s(ista) = var_fcst_3d(i_g,j_g,k_g)
!                       endif
 
                        sumobs = sumobs + var_s(ista)
                        sumanl = sumanl + var_fcst_s(ista)
                        cnt = cnt + 1.

                      else
                        iwrite = iwrite + 1
                        if(iwrite .le. 50)then
                          write(6,*)' WARNING: k is outside domain '
     1                             ,ista,k_g       
                        endif

                      endif

1112                endif ! ob is in domain
                  endif ! non-missing observation (valid value)
                enddo ! ista

                write(6,*)
                write(6,*)' Generic stats, cnt = ',nint(cnt)

                if(nint(cnt) .eq. 0)then
                  write(6,*)' Skipping call to stats_1d'
                  xbar = rmiss
                  ybar = rmiss
                  std = rmiss
                  goto 900
                endif

1200            if(var_2d .eq. 'U3')then
                  call stats_1d(maxsta,var_fcst_s,var_s 
     1                   ,'U Wind Component:       '
     1                   ,a_t,b_t,xbar,ybar
     1                   ,bias,std,r_missing_data,istatus)
                  if(istatus .eq. 1)then
                      rms_a(ifield,imodel,itime_fcst) = std
                  else
                      rms_a(ifield,imodel,itime_fcst) = r_missing_data
                  endif
                elseif(var_2d .eq. 'V3')then
                  call stats_1d(maxsta,var_fcst_s,var_s
     1                   ,'V Wind Component:       '
     1                   ,a_t,b_t,xbar,ybar
     1                   ,bias,std,r_missing_data,istatus)
                  if(istatus .eq. 1)then
                      rms_a(ifield,imodel,itime_fcst) = std
                  else
                      rms_a(ifield,imodel,itime_fcst) = r_missing_data
                  endif
                elseif(var_2d .eq. 'W3')then
                  xbar = rms_a(1,imodel,itime_fcst)
                  ybar = rms_a(2,imodel,itime_fcst)
                  if(xbar .ne. r_missing_data .AND.
     1               ybar .ne. r_missing_data       )then
                      std = sqrt(rms_a(1,imodel,itime_fcst)**2
     1                          +rms_a(2,imodel,itime_fcst)**2)
                  else
                      xbar = rmiss
                      ybar = rmiss
                      std = rmiss                         
                  endif
                  rms_a(ifield,imodel,itime_fcst) = std
                elseif(var_2d .eq. 'T3')then
                  call stats_1d(max_obs,var_fcst_s,var_s
     1                   ,'Temperature:            '
     1                   ,a_t,b_t,xbar_k,ybar_k
     1                   ,bias,std,r_missing_data,istatus)
                  if(istatus .eq. 1)then
                      xbar = k_to_c(xbar_k)
                      ybar = k_to_c(ybar_k)
                  else
                      xbar = rmiss
                      ybar = rmiss
                      std = rmiss                         
                  endif
                  rms_a(ifield,imodel,itime_fcst) = std
                endif

900             continue

990             continue

                write(6,*)
                write(6,*)' Writing to lun_out ',lun_out
                call cv_i4tim_asc_lp(i4_valid,atime_s,istatus)
                write(6,710)atime_s,xbar,ybar,std
                write(lun_out,710)atime_s,xbar,ybar,std
710             format(1x,a24,3f10.3)

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
