          
        subroutine verif_fcst_grid_2d(i4time_sys,a9time,laps_cycle_time,       
     1                  ni,nj,
     1                  nk,
     1                  maxsta,
     1                  r_missing_data,
     1                  model_verif_intvl,
     1                  n_fcst_times,
     1                  j_status)

        use mem_namelist, ONLY: path_to_gps

        real var_anal_2d(ni,nj)
        real var_fcst_2d(ni,nj)
        real var_prst_2d(ni,nj)
        real u_fcst_2d(ni,nj)
        real v_fcst_2d(ni,nj)
        real lat(ni,nj)
        real lon(ni,nj)
        real topo(ni,nj)
        real rlaps_land_frac(ni,nj)

        real k_to_f

        logical lmask_and_3d(ni,nj,nk)
        logical lmask_rqc_3d(ni,nj,nk)
        logical l_parse, l_parse_result

        integer       maxbgmodels
        parameter     (maxbgmodels=10)
        character*30  c_fdda_mdl_src(maxbgmodels)

        integer max_fcst_times
        parameter (max_fcst_times=200)

        character EXT*31, directory*255, c_model*30

        character*10  units_2d
        character*125 comment_2d
        character*3 var_2d,var_2d_wind
        character*9 a9time,a9time_valid,a9time_init
        character*24 atime_s
        character*150 hist_dir, cont_dir, verif_dir
        character*150 hist_file

        integer n_fields
        parameter (n_fields=1)
        character*10 ext_anal_a(n_fields), ext_fcst_a(n_fields)
        character*10 var_a(n_fields)

        logical l_persist, l_good_persist 

!       Specify what is being verified
        data ext_fcst_a /'fsf'/ ! 2-D
        data ext_anal_a /'lcv'/ ! 2-D
        data var_a /'S8A'/
      
        real rms_a (n_fields,maxbgmodels,0:max_fcst_times)

!       Declarations from compare_analysis_to_rad
        real cloud_frac_vis_a(ni,nj),tb8_k(ni,nj),t_gnd_k(ni,nj)
     1        ,t_sfc_k(ni,nj),cvr_max(ni,nj),cvr_sao_max(ni,nj)
     1        ,dbz_max_2d(ni,nj),solar_alt(ni,nj),swi_2d(ni,nj)

!       How much the solar radiation varies with changes in cloud fraction
        real cvr_scl_a(ni,nj) 

        real rad_clr(ni,nj)

        real airmass_2d(ni,nj)
        real pw_2d(ni,nj)
        real trans_h2o_2d(ni,nj)

        real dum_2d(ni,nj)

        character c_stations(maxsta)*3
        character stn(maxsta)*20

!       real lat_s(maxsta), lon_s(maxsta), elev_s(maxsta)
        real var_s(maxsta)
        real rad2_s(maxsta)
        real swi_s(maxsta)
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

        write(6,*)' Start subroutine verif_fcst_grid_2d...'

        l_persist = .true. ! Add persistence forecast for evaluation

        rmiss = -99.9

        call get_laps_domain_95(ni,nj,lat,lon,topo
     1           ,rlaps_land_frac,grid_spacing_cen_m,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error getting LAPS domain'
            return
        endif

!       n_fcst_times = 2 ! 38

        rms_a = -99.9

        thresh_var = 20. ! lowest threshold for this variable

        i4_initial = i4time_sys

        lun_in = 21

!       Get fdda_model_source from static file
        call get_fdda_model_source(c_fdda_mdl_src,n_fdda_models,istatus)

        if(l_persist .eqv. .true.)then
            n_fdda_models = n_fdda_models + 1
            c_fdda_mdl_src(n_fdda_models) = 'persistence'
            write(6,*)' Adding persistence to fdda_models'
        endif

        write(6,*)' n_fdda_models = ',n_fdda_models
        write(6,*)' c_fdda_mdl_src = '
     1            ,(c_fdda_mdl_src(m),m=1,n_fdda_models)

!       if(n_fdda_models .ne. n_models + 1)then
!           write(6,*)' ERROR n_models differs from n_fdda_models '
!    1                       ,n_models,n_fdda_models
!           stop
!       endif

        if(n_fdda_models .lt. 2)then
            write(6,*)' WARNING: n_fdda_models is less than 2'
            write(6,*)' Check nest7grid.parms specification'
        endif

        do ifield = 1,n_fields

          var_prst_2d = r_missing_data
          l_good_persist = .false.

          var_2d = var_a(ifield)(1:3)
          call s_len(var_2d,lenvar)

          write(6,*) 
          write(6,*) 
          write(6,*)' Processing field ',trim(var_2d) 

          call get_directory('verif',verif_dir,len_verif)

          if(trim(var_2d) .eq. 'S8A')then      
              istart_fcst = 1
          else
              istart_fcst = 0
          endif

          lun_out = 39

          hist_dir = verif_dir(1:len_verif)//var_2d(1:lenvar)
     1                                     //'/pt'
!    1                                     //c_model(1:len_model)
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

            call s_len(c_model,len_model)

            cont_dir = verif_dir(1:len_verif)//var_2d(1:lenvar)
     1                                       //'/cont/'
     1                                       //c_model(1:len_model)
     1                                       //'/'
            len_cont = len_verif + 6 + lenvar + len_model

            if(c_model(1:3) .ne. 'lga')then

              write(6,*)
              write(6,*)' Processing model ',c_model

              call s_len(c_model,len_model)

!             if(istart_fcst .eq. 1)then
!                 call cv_i4tim_asc_lp(i4_initial,atime_s,istatus)
!                 write(lun_out,710)atime_s,rmiss,rmiss,rmiss
!             endif

              do itime_fcst = 0,n_fcst_times

                cnt = 0.

                itime = itime_fcst + 1

                i4_valid = i4_initial 
     1                   + itime_fcst * model_verif_intvl

                call make_fnam_lp(i4_valid,a9time_valid,istatus)

                write(6,*) 
                write(6,*)' Processing time ',a9time_valid

                call cv_i4tim_asc_lp(i4_valid,atime_s,istatus)

                xbar = rmiss
                ybar = rmiss
                std = rmiss

                write(6,*)' Reading 2D analysis i4_valid = ',i4_valid

                i4_fcst = i4_valid - i4_initial

                ext = ext_anal_a(ifield)

                call get_laps_2d(i4_valid,ext,var_2d,units_2d       
     1              ,comment_2d,ni,nj,var_anal_2d,istatus)       

                if(istatus .eq. -1)then
                    nmissing = 0
                    do i = 1,ni
                    do j = 1,nj
                        if(var_anal_2d(i,j) .eq. r_missing_data)then
                            nmissing = nmissing + 1
                        endif
                    enddo ! j
                    enddo ! i

                    pct_missing = (float(nmissing) / float(ni*nj))*100.

                    if(pct_missing .le. 10.)then
                        write(6,*)' WARNING: missing analysis data for '
     1                             ,var_2d,nmissing,pct_missing
                    else
                        write(6,*)' ERROR: missing analysis data for '
     1                             ,var_2d,nmissing,pct_missing
                        goto 980
                    endif

                elseif(istatus .ne. 1)then
                    write(6,*)' Error reading 2D Analysis for '
     1                         ,var_2d
                    goto 980
                endif

                write(6,*)var_2d,' anal range is ',minval(var_anal_2d)
     1                                            ,maxval(var_anal_2d)

                if(itime_fcst .eq. 0)then
                    var_prst_2d = var_anal_2d
                    l_good_persist = .true.
                endif

                if(c_fdda_mdl_src(imodel) .ne. 'persistence')then

                 if(itime_fcst .lt. istart_fcst)then
                  write(6,*)' Forecast unavailable at time ',itime_fcst       
                  goto 990
                 endif

                 if(trim(var_2d) .ne. 'SSF')then

                  write(6,*)' Reading forecast field ',atime_s

!                 Read forecast field
                  ext = ext_fcst_a(ifield)
                  call get_directory(ext,directory,len_dir)

                  DIRECTORY=directory(1:len_dir)//c_model(1:len_model)
     1                                          //'/'

                  call get_lapsdata_2d(i4_initial,i4_valid
     1                          ,directory,var_2d
     1                          ,units_2d,comment_2d
     1                          ,ni,nj
     1                          ,var_fcst_2d
     1                          ,istatus)

!                 Suppress 00hr SWI if from some models since it would have been
!                 pulled in via LFMPOST from a LAPS analysis

                  if(trim(var_2d) .eq. 'SWI')then              
                     l_parse_result=
     1                   (l_parse(c_model(1:len_model),'hrrr') .OR.
     1                    l_parse(c_model(1:len_model),'rr'  ) .OR.
     1                    l_parse(c_model(1:len_model),'rap-nh'))
                     if(l_parse_result .eqv. .true.)then       
                        if(i4_valid .eq. i4_initial)then
                           write(6,*)' Suppressing 00hr SWI from '
     1                                ,c_model(1:len_model)
                          istatus = 0
                        endif
                     endif
                  endif

                  if(istatus .ne. 1)then
                       write(6,*)' Error reading 2D Forecast for '
     1                           ,var_2d
                       goto 990
                  endif

                  write(6,*)var_2d,' fcst range is ',minval(var_fcst_2d)
     1                                              ,maxval(var_fcst_2d)

                 else ! SSF

                  write(6,*)' Reading forecast U/V to yield speed'

!                 Read forecast field
                  ext = ext_fcst_a(ifield)
                  call get_directory(ext,directory,len_dir)

                  DIRECTORY=directory(1:len_dir)//c_model(1:len_model)
     1                                          //'/'

                  var_2d_wind = 'USF'
                  call get_lapsdata_2d(i4_initial,i4_valid
     1                          ,directory,var_2d_wind
     1                          ,units_2d,comment_2d
     1                          ,ni,nj
     1                          ,u_fcst_2d
     1                          ,istatus)
                  if(istatus .ne. 1)then
                       write(6,*)' Error reading 2D Forecast for '
     1                           ,var_2d_wind
                       goto 990
                  endif

                  var_2d_wind = 'VSF'
                  call get_lapsdata_2d(i4_initial,i4_valid
     1                          ,directory,var_2d_wind
     1                          ,units_2d,comment_2d
     1                          ,ni,nj
     1                          ,u_fcst_2d
     1                          ,istatus)
                  if(istatus .ne. 1)then
                       write(6,*)' Error reading 2D Forecast for '
     1                           ,var_2d_wind
                       goto 990
                  endif

                  var_2d_wind = 'VSF'
                  call get_lapsdata_2d(i4_initial,i4_valid
     1                          ,directory,var_2d_wind
     1                          ,units_2d,comment_2d
     1                          ,ni,nj
     1                          ,v_fcst_2d
     1                          ,istatus)
                  if(istatus .ne. 1)then
                       write(6,*)' Error reading 2D Forecast for '
     1                           ,var_2d_wind
                       goto 990
                  endif

!                 Convert U and V to wind speed
                  do i = 1,ni
                  do j = 1,nj
                      call uv_to_disp(u_fcst_2d(i,j)
     1                               ,v_fcst_2d(i,j)
     1                               ,dir_dum                 
     1                               ,var_fcst_2d(i,j))
                  enddo ! j
                  enddo ! i

                 endif ! SSF field

                elseif(l_good_persist .eqv. .true.)then
                    write(6,*)
     1                  ' Setting forecast to persistence gridded data '
     1                  ,var_2d
                    var_fcst_2d = var_prst_2d
                else
                    write(6,*)' Persistence fcst unavailable'
                    goto 990
                endif

                if(trim(var_2d) .eq. 'SWI')then
                  swi_2d = var_fcst_2d
                endif
               
                ista = 0

                do i = 1,ni
                do j = 1,nj
                  ista = ista + 1

!                 write(6,*)ista,var_s(ista)
                  swi_s(ista) = r_missing_data 
                  var_fcst_s(ista) = r_missing_data 
                  rad2_s(ista) = r_missing_data 

                  var_s(ista) = var_anal_2d(i,j)

                  if(var_s(ista) .gt. threshval .AND.
     1               var_s(ista) .ne. r_missing_data)then ! valid ob  
                
                    i_i = i
                    i_j = j

                    ii_s(ista) = i_i
                    jj_s(ista) = i_j

                    if(.true.)then

                      var_fcst_s(ista) = var_fcst_2d(i_i,i_j)            

                      if(cnt .le. 50.)then
                        write(6,1001)var_2d,ista,i_i,i_j 
     1                              ,var_s(ista),var_fcst_s(ista)
1001                    format(1x,a3,' ob/fcst ',3i8,2f9.3)
                      endif

                      sumobs = sumobs + var_s(ista)
                      sumanl = sumanl + swi_2d(i_i,i_j)
                      sumcld = sumcld + cvr_max(i_i,i_j)
                      sumalt = sumalt + solar_alt(i_i,i_j)
                      sumscl = sumscl + cvr_scl_a(i_i,i_j)
                      cnt = cnt + 1.

                      cvr_s(ista) = cvr_max(i_i,i_j)
1112                  continue

                    endif 

                  else
                    var_fcst_s(ista) = r_missing_data

                  endif ! valid ob                                  

                enddo ! j
                enddo ! i

                write(6,*)
                write(6,*)' Generic stats, cnt = ',nint(cnt)

                if(cnt .eq. 0.)then
                  write(6,*)' No obs - skipping calculation of stats'
                  goto 980 ! Skip calculation of stats and write initialized flag values
                endif
 
1200            if(trim(var_2d) .eq. 'S8A')then
                  maxsta = ni*nj
                  call stats_1d(maxsta,var_fcst_s,var_s
     1                   ,'11 micron Sat: '
     1                   ,a_t,b_t,xbar,ybar
     1                   ,bias,std,r_missing_data,istatus)
                  rms_a(ifield,imodel,itime_fcst) = std
                endif

980             continue

990             continue
 
!               In General:
!               xbar is mean forecast value at the stations
!               ybar is mean observed value at the stations

!               For WSF (Surface Wind)
!               xbar is the rms of U
!               ybar is the rms of V

                write(6,*)
                write(6,*)' Writing to lun_out ',lun_out
                write(6,710)atime_s,xbar,ybar,std
                write(lun_out,710)atime_s,xbar,ybar,std
!               write(39,710)atime_s,xbar,ybar,std
!               write(38,710)atime_s,xbar,ybar,std
710             format(1x,a24,3f10.3)

              enddo ! itime_fcst

              write(lun_out,*)
              write(lun_out,*)

            endif ! c_model .ne. lga

          enddo ! model

          write(6,*)' Closing lun_out ',lun_out
          close (lun_out) 

        enddo ! fields

 999    write(6,*)' End of subroutine verif_fcst_grid_2d'

        return

        end
