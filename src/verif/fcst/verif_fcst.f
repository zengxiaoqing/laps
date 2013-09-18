  
        program verif_fcst_main

        use mem_namelist, ONLY: read_namelist_laps
        use mem_namelist, ONLY: model_fcst_intvl,model_cycle_time
     1                                          ,model_fcst_len

        character*150 static_dir,filename
        character*9 a9time
	character*300 dir_t,filenamet
        character*10 c_n_fcst_times,c_model_fcst_intvl
        character*150 verif_dir, n_plot_times_file
        integer, parameter :: lun=120

!       This can be changed to read the modeltime for the forecast
!       initialization?
!       call get_systime(i4time,a9time,istatus)
!       if(istatus .ne. 1)go to 999

        ISTAT = init_timer()
        I4_elapsed = ishow_timer()

        call get_directory('time',dir_t,istatus)
        call s_len(dir_t,len_dir_t)
        filenamet = dir_t(1:len_dir_t)//'/modelvtime.dat'
        write(6,*)' filenamet = ',trim(filenamet)
        open(lun,file=filenamet,status='old')
        read(lun,*)a9time
        close(lun)
        call i4time_fname_lp(a9time,i4time,istatus) 
        write(6,*)' model time = ',a9time

        call get_grid_dim_xy(NX_L,NY_L,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'Error getting horizontal domain dimensions'
           go to 999
        endif

        call get_laps_dimensions(NZ_L,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'Error getting vertical domain dimension'
           go to 999
        endif

        call get_r_missing_data(r_missing_data,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'Error getting r_missing_data'
           go to 999
        endif

        call get_laps_cycle_time(laps_cycle_time,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'Error getting laps_cycle_time'
           go to 999
        endif

!       Read moisture parameters into module memory structure
        call get_directory('static',static_dir,len_dir)
        filename = static_dir(1:len_dir)//'/moisture_switch.nl'
        call read_namelist_laps('moisture_anal',filename)

        maxsta = 1000000
        max_obs = 1000000

!       Use getarg for model forecast interval
!       call getarg(1,c_model_fcst_intvl)
!       read(c_model_fcst_intvl,*)model_fcst_intvl

!       model_fcst_intvl = 3600  
!       model_fcst_intvl = 10800

!       Check getarg for number of forecasts (superseded by parms calculation)
!       call getarg(2,c_n_fcst_times)
!       read(c_n_fcst_times,*)n_fcst_times

        model_verif_intvl = max(laps_cycle_time,model_fcst_intvl)
        n_fcst_times = (model_fcst_len*60) / model_verif_intvl

        write(6,*)
     1  ' model_fcst_intvl (via namelist) / n_fcst_times (from parms) ='      
     1             ,model_fcst_intvl,n_fcst_times

        if(.true.)then

            write(6,*)' Calling verif_fcst_pt_2d'
          
            call verif_fcst_pt_2d(i4time,a9time,laps_cycle_time,
     1                     NX_L,NY_L,
     1                     NZ_L,
     1                     maxsta,
     1                     r_missing_data,
     1                     model_verif_intvl,
     1                     n_fcst_times,
     1                     j_status)

            I4_elapsed = ishow_timer()

            write(6,*)' Calling verif_fcst_pt_3d'

            call verif_fcst_pt_3d(i4time,a9time,laps_cycle_time,
     1                     NX_L,NY_L,
     1                     NZ_L,
     1                     maxsta,max_obs,
     1                     r_missing_data,
     1                     model_verif_intvl,
     1                     n_fcst_times,
     1                     j_status)

            I4_elapsed = ishow_timer()

        endif

        maxsta=NX_L*NY_L ! some 2D grids are fed into 1-D arrays

        call verif_fcst_grid_2d(i4time,a9time,laps_cycle_time,
     1                     NX_L,NY_L,
     1                     NZ_L,
     1                     maxsta,
     1                     r_missing_data,
     1                     model_verif_intvl,
     1                     n_fcst_times,
     1                     j_status)

        I4_elapsed = ishow_timer()

!       Read n_plot_times from file
        call get_directory('verif',verif_dir,len_verif)
        lun_plot_times = 42
        n_plot_times_file = verif_dir(1:len_verif)//'/n_fcst_times.dat'
        open(lun_plot_times,file=n_plot_times_file,status='old')
        read(lun_plot_times,*)n_plot_times
        close(lun_plot_times)

        call verif_fcst_composite(i4time,a9time,                    
     1                     model_fcst_intvl,
     1                     model_fcst_len,
     1                     model_cycle_time,
     1                     laps_cycle_time,
     1                     r_missing_data,
!    1                     model_verif_intvl,
     1                     n_plot_times,
     1                     j_status)

        I4_elapsed = ishow_timer()

999     continue

        end
          
        subroutine verif_fcst_pt_2d(i4time_sys,a9time,laps_cycle_time,
     1                  ni,nj,
     1                  nk,
     1                  maxsta,
     1                  r_missing_data,
     1                  model_verif_intvl,
     1                  n_fcst_times,
     1                  j_status)

        use mem_namelist, ONLY: path_to_gps, iverbose

        include 'read_sfc.inc'

        real var_anal_2d(ni,nj)
        real var_fcst_2d(ni,nj)
        real var_fcst_2d_prev(ni,nj)
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
        logical l_parse, l_parse_result, l_exist

        integer       maxbgmodels
        parameter     (maxbgmodels=10)
        character*30  c_fdda_mdl_src(maxbgmodels)

        integer max_fcst_times
        parameter (max_fcst_times=200)

        integer max_regions
        parameter (max_regions=10)

        integer il(maxbgmodels,max_fcst_times,max_regions)
        integer ih(maxbgmodels,max_fcst_times,max_regions)
        integer jl(maxbgmodels,max_fcst_times,max_regions)
        integer jh(maxbgmodels,max_fcst_times,max_regions)

        character EXT*31, directory*255, c_model*30

        character*10  units_2d
        character*125 comment_2d
        character*3 var_2d,var_2d_rto,var_2d_wind
        character*9 a9time,a9time_valid,a9time_init
!       character*24 atime_s
        character*150 hist_dir, cont_dir, verif_dir
        character*150 hist_file

        integer n_fields
        parameter (n_fields=13)
        character*10 ext_anal_a(n_fields), ext_fcst_a(n_fields)
        character*10 var_a(n_fields)
        integer ipersist_a(n_fields) ! persistence flag for each field              
        logical l_persist, l_good_persist

!       Specify what is being verified
        data ext_fcst_a 
     1     /'fsf','fsf','fsf','fsf','fsf','fsf','   ','fsf','fsf','fsf',
     1      'fsf','fsf','fsf'/ ! 2-D
!       data ext_anal_a /'lps'/ ! 3-D reflectivity
        data var_a      
     1     /'SWI','TSF','DSF','USF','VSF','SSF','WSF','TPW','R01','RTO',
     1      'R03','R06','R24'/
      
        data ipersist_a /0,0,0,0,0,0,0,1,0,0,0,0,0/        

        real rms_a (n_fields,maxbgmodels,0:max_fcst_times)
        real cnt_a (n_fields,maxbgmodels,0:max_fcst_times)

        integer contable(0:1,0:1)

        integer maxthr
        parameter (maxthr=5)

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

!       GPS obs
        integer gps_n, gps_indomain
        parameter (gps_n = 1000000)
        real gps_tpw(gps_n)
        real gps_wet(gps_n)
        real gps_error(gps_n)
        real gps_xy(2,gps_n)
        real gps_elv(gps_n)
        real gps_tim(gps_n)
!       character*256 path_to_gps

!       End Declarations 

        write(6,*)' Start subroutine verif_fcst_pt_2d...'

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

        do ifield = 1,n_fields

          if(ipersist_a(ifield) .eq. 1)then
            l_persist = .true.
          else
            l_persist = .false.
          endif

          var_prst_2d = r_missing_data
          l_good_persist = .false.

!         Get fdda_model_source from static file
          call get_fdda_model_source(c_fdda_mdl_src,n_fdda_models
     1                                             ,istatus)

          if(l_persist .eqv. .true.)then
              n_fdda_models = n_fdda_models + 1
              c_fdda_mdl_src(n_fdda_models) = 'persistence'
              write(6,*)' Adding persistence to fdda_models'
          endif

          write(6,*)' n_fdda_models = ',n_fdda_models
          write(6,*)' c_fdda_mdl_src = '
     1              ,(c_fdda_mdl_src(m),m=1,n_fdda_models)

!         if(n_fdda_models .ne. n_models + 1)then
!           write(6,*)' ERROR n_models differs from n_fdda_models '
!    1                       ,n_models,n_fdda_models
!           stop
!         endif

          if(n_fdda_models .lt. 2)then
              write(6,*)' WARNING: n_fdda_models is less than 2'
              write(6,*)' Check nest7grid.parms specification'
          endif

          var_2d = var_a(ifield)(1:3)
          call s_len(var_2d,lenvar)

          write(6,*) 
          write(6,*) 
          write(6,*)' Processing field ',trim(var_2d) 

          call get_directory('verif',verif_dir,len_verif)

          if(var_2d(1:2) .eq. 'R0' .OR. var_2d(1:2) .eq. 'R2')then      
              istart = 1
          else
              istart = 0
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

              if(istart .eq. 1)then
                  call cv_i4tim_asc_lp(i4_initial,atime_s,istatus)
                  cnt = 0.
                  write(lun_out,710)atime_s,rmiss,rmiss,rmiss,nint(cnt)
              endif

              do itime_fcst = istart,n_fcst_times

                cnt = 0.
                cnt_write = 0.

                itime = itime_fcst + 1

                i4_valid = i4_initial 
     1                   + itime_fcst * model_verif_intvl

                call make_fnam_lp(i4_valid,a9time_valid,istatus)

                write(6,*) 
                write(6,*)' Processing time ',a9time_valid

                if(trim(var_2d) .eq. 'WSF')then
                    call cv_i4tim_asc_lp(i4_valid,atime_s,istatus)
                    goto 1200
                endif

                xbar = rmiss
                ybar = rmiss
                std = rmiss

                write(6,*)' Reading surface obs - i4_valid = ',i4_valid

                i4_fcst = i4_valid - i4_initial

                if(trim(var_2d) .ne. 'TPW')then
!                 Read surface obs (regular LSO)
                  call read_surface_data(
     &         i4_valid,atime_s,n_obs_g,n_obs_b,
     &         obstime,wmoid,stations,provider,wx_s,reptype,autostntype,
     &         lat_s,lon_s,elev_s,t_s,td_s,rh_s,dd_s,ff_s,ddg_s,ffg_s,
     &         alt_s,pstn_s,pmsl_s,delpch,delp,vis_s,solar_s,sfct,sfcm,
     &         pcp1,pcp3,pcp6,pcp24,snow,kloud_s,max24t,min24t,t_ea,
     &         td_ea,rh_ea,dd_ea,ff_ea,alt_ea,p_ea,vis_ea,solar_ea,
     &         sfct_ea,sfcm_ea,pcp_ea,snow_ea,store_amt,store_hgt,
     &         maxsta,jstatus)

                  if(jstatus .ne. 1)then
                      call cv_i4tim_asc_lp(i4_valid,atime_s,istatus)
                      write(6,*)' Could not read surface obs at: ',
     1                          atime_s
                      goto 990
                  endif

                  write(6,*)' Number of obs in box  ',n_obs_b,atime_s

                else
                  call cv_i4tim_asc_lp(i4_valid,atime_s,istatus)

                  write(6,*)' call read_gps_obs: ',trim(path_to_gps)
                  inquire(file=trim(path_to_gps),exist=l_exist)
                  if(l_exist .eqv. .true.)then
                      write(6,*)' Path to GPS from namelist exists'         
                  else
                      write(6,*)
     1                  ' Path to GPS non-existent, trying /public'
                      path_to_gps = '/public/data/gpsmet/netcdf/'
                  endif
                  lun_hmg = 0
                  i4beg = i4_valid - 960
                  i4end = i4_valid + 840
                  istatus = 0
                  bad_sfc = 0.0
                  gps_missing = -9.99
                  call read_gps_obs (lun_hmg, path_to_gps, i4beg, i4end,
     1                    ni, nj, lat, lon, bad_sfc,
     1                    gps_tpw, gps_wet, gps_error, gps_xy, gps_elv,
     1                    gps_tim, gps_indomain, gps_n, istatus)
                  if(istatus .eq. 1)then
                    write(6,*)
     1                    ' Success reading GPS obs, #obs in domain = '
     1                     ,gps_indomain
                    if(gps_indomain .eq. 0.)then
                      write(6,*)' No GPS obs - skip stats calculation'
                      goto 980 ! Skip calculation of stats and write initialized flag values
                    endif
                  else
                    write(6,*)' Failure reading GPS obs at:',atime_s
                    goto 980
                  endif

                endif

                if(trim(var_2d) .eq. 'SWI')then
                  clear_sky_ghi = r_missing_data
                  do i = 1,maxsta  
                      call solalt(lat_s(i),lon_s(i),i4_valid,sol_alt)
                      call qc_solar_ob(solar_s(i),sol_alt 
     1                                ,r_missing_data,iqc,clear_sky_ghi)

                      call get_sfc_obtime(obstime(i),i4_valid,i4time_ob        
     1                                   ,istatus)
                      if(abs(i4time_ob - i4_valid) .gt. 60)then
                          iqc = 1 ! ob is more than a minute from valid time
                      endif

                      if(iqc .eq. 0)then
                          var_s(i) = solar_s(i)
                      else
                          var_s(i) = r_missing_data
                          if(iverbose .eq. 1)then
                              write(6,*)' Solar ob QCd out '
     1                        ,i,stations(i),sol_alt,solar_s(i)
     1                        ,provider(i),obstime(i)       
                          endif
                      endif
                  enddo ! i
                  threshval = 0.
                elseif(trim(var_2d) .eq. 'TSF')then
                  var_s = t_s
                  threshval = -99.9
                elseif(trim(var_2d) .eq. 'DSF')then
                  var_s = td_s
                  threshval = -99.9
                elseif(trim(var_2d) .eq. 'USF')then
                  threshval = -99.9
                  do i = 1,maxsta  
                    if(dd_s(i) .ne. threshval .and. 
     1                 ff_s(i) .ne. threshval)then
                      call disp_to_uv(dd_s(i),ff_s(i),u_s,v_s)
                      var_s(i) = u_s * 0.518 ! convert to m/s
                    else
                      var_s(i) = threshval
                    endif
                  enddo ! i
                elseif(trim(var_2d) .eq. 'VSF')then
                  threshval = -99.9
                  do i = 1,maxsta  
                    if(dd_s(i) .ne. threshval .and. 
     1                 ff_s(i) .ne. threshval)then
                      call disp_to_uv(dd_s(i),ff_s(i),u_s,v_s)
                      var_s(i) = v_s * 0.518 ! convert to m/s
                    else
                      var_s(i) = threshval
                    endif
                  enddo ! i
                elseif(trim(var_2d) .eq. 'SSF')then
                  threshval = -99.9
                  do i = 1,maxsta
                    if(ff_s(i) .ne. threshval)then
                      var_s(i) = ff_s(i) * 0.518 ! convert to m/s
                    else
                      var_s(i) = threshval
                    endif
                  enddo ! i
                elseif(trim(var_2d) .eq. 'TPW')then
                  threshval = gps_missing
                  do i = 1,maxsta
                    if(i .le. gps_indomain)then
                      var_s(i) = gps_tpw(i)
                    else
                      var_s(i) = r_missing_data
                    endif
                  enddo ! i
                elseif(trim(var_2d) .eq. 'R01_old')then
                  threshval = -99.9
                  if(model_verif_intvl .eq. 3600)then
                      var_s = pcp1
                  elseif(model_verif_intvl .eq. 10800)then
                      var_s = pcp3
                  elseif(model_verif_intvl .eq. 21600)then
                      var_s = pcp6
                  elseif(model_verif_intvl .eq. 86400)then
                      var_s = pcp24
                  else
                      var_s = r_missing_data  
                  endif
                elseif(trim(var_2d) .eq. 'RTO')then
                  threshval = -99.9
                  if(i4_fcst .eq. 3600)then
                      var_s = pcp1
                  elseif(i4_fcst .eq. 10800)then
                      var_s = pcp3
                  elseif(i4_fcst .eq. 21600)then
                      var_s = pcp6
                  elseif(i4_fcst .eq. 86400)then
                      var_s = pcp24
                  else
                      var_s = r_missing_data  
                  endif
                elseif(trim(var_2d) .eq. 'R01')then
                  threshval = -99.9
                  var_s = pcp1
                elseif(trim(var_2d) .eq. 'R03')then
                  threshval = -99.9
                  var_s = pcp3
                elseif(trim(var_2d) .eq. 'R06')then
                  threshval = -99.9
                  var_s = pcp6
                elseif(trim(var_2d) .eq. 'R24')then
                  threshval = -99.9
                  var_s = pcp24
                endif

!               Read forecast field
                if(c_fdda_mdl_src(imodel) .ne. 'persistence')then
                 if(trim(var_2d) .eq. 'R01'  .OR.
     1              trim(var_2d) .eq. 'R03'  .OR.    
     1              trim(var_2d) .eq. 'R06'  .OR.    
     1              trim(var_2d) .eq. 'R24'       )then

                    if(trim(var_2d) .eq. 'R01')i4_acc =  3600
                    if(trim(var_2d) .eq. 'R03')i4_acc = 10800
                    if(trim(var_2d) .eq. 'R06')i4_acc = 21600
                    if(trim(var_2d) .eq. 'R24')i4_acc = 86400

                    if(i4_initial .eq. (i4_initial/i4_acc)*i4_acc)then
                      if(i4_valid .eq. (i4_valid/i4_acc)*i4_acc)then
                        var_2d_rto = 'RTO'
                        call get_directory(ext,directory,len_dir)
                        DIRECTORY=
     1                  directory(1:len_dir)//c_model(1:len_model)//'/'

                        call get_lapsdata_2d(i4_initial,i4_valid
     1                          ,directory,var_2d_rto
     1                          ,units_2d,comment_2d
     1                          ,ni,nj
     1                          ,var_fcst_2d
     1                          ,istatus)

                        if(istatus .ne. 1)then
                          write(6,*)' Error reading 2D Forecast for '
     1                           ,var_2d
                          goto 990
                        endif

                        call get_lapsdata_2d(i4_initial,i4_valid-i4_acc
     1                          ,directory,var_2d_rto
     1                          ,units_2d,comment_2d
     1                          ,ni,nj
     1                          ,var_fcst_2d_prev
     1                          ,istatus)

                        if(istatus .ne. 1)then
                          write(6,*)
     1                    ' Error reading 2D Prev Forecast for ',var_2d
                          goto 990
                        endif

                        write(6,*)
     1                    ' Subtracting precip totals over interval'
                        var_fcst_2d = var_fcst_2d - var_fcst_2d_prev

                      else
                        goto 990
                      endif ! useable valid time

                    else
                      goto 990
                    endif ! useable initial time

                 elseif(trim(var_2d) .ne. 'SSF')then

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

                  if(istatus .ne. 1)then
                       write(6,*)' Error reading 2D Forecast for '
     1                           ,var_2d
                       goto 990
                  endif

                  write(6,*)var_2d,' range is ',minval(var_fcst_2d)
     1                                         ,maxval(var_fcst_2d)

                  if(trim(var_2d) .eq. 'TPW')then
                      if(maxval(var_fcst_2d) .eq. 0.)then
                          write(6,*)' QC issue with TPW forecast field' 
                          goto 990
                      endif
                  endif

                  if(l_persist .eqv. .true. .and. 
     1               l_good_persist .eqv. .false. .and.
     1               itime_fcst .eq. 0                  )then
                       write(6,*)imodel,l_persist,l_good_persist
     1                          ,itime_fcst,istatus
     1                          ,' Setting persistence to 00 hr fcst '
     1                          ,var_2d
                       var_prst_2d = var_fcst_2d

!                      QC check for TPW persistence forecast
                       if(maxval(var_prst_2d) .eq. 0. .and.
     1                    minval(var_prst_2d) .eq. 0. .and.
     1                    trim(var_2d) .eq. 'TPW'     )then
                           write(6,*)
     1                    ' WARNING, persistence appears to be all zero'
                       else
                           l_good_persist = .true.
                       endif

!                 else
!                      write(6,*)imodel,l_persist,l_good_persist
!    1                          ,itime_fcst,istatus,
!    1                         ' Not setting persistence to 00 hr fcst'
                  endif

!                 Suppress 00hr SWI if from wrf-hrrr since it would have been
!                 pulled in via LFMPOST from a LAPS analysis

                  if(trim(var_2d) .eq. 'SWI')then              
                     l_parse_result=l_parse(c_model(1:len_model),'hrrr')
                     if(l_parse_result .eqv. .true.)then       
                        if(i4_valid .eq. i4_initial)then
                           write(6,*)' Suppressing 00hr hrrr SWI'
                          istatus = 0
                        endif
                     endif
                  endif

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

                 endif ! Variable of interest (precip,ssf,generic)

                elseif(l_good_persist .eqv. .true.)then
                    write(6,*)' Setting forecast to persistence ',var_2d
                    var_fcst_2d = var_prst_2d
                else
                    write(6,*)' Persistence fcst unavailable'
                    goto 990
                endif

                if(trim(var_2d) .eq. 'SWI')then
                  swi_2d = var_fcst_2d
                endif

                do ista = 1,maxsta  
!                 stn(ista) = c_stations(ista)(1:3)
                  stn(ista) = stations(ista)(1:3)

!                 write(6,*)ista,var_s(ista)
                  swi_s(ista) = r_missing_data 
                  var_fcst_s(ista) = r_missing_data 
                  rad2_s(ista) = r_missing_data 

                  if(var_s(ista) .gt. threshval .AND.
     1               var_s(ista) .ne. r_missing_data)then ! valid ob  
                
                    if(trim(var_2d) .ne. 'TPW')then
                      call latlon_to_rlapsgrid(lat_s(ista),lon_s(ista)
     1                          ,lat,lon,ni,nj,ri,rj,istatus)
                    else
                      ri = gps_xy(1,ista)
                      rj = gps_xy(2,ista)
                    endif

                    i_i = nint(ri)
                    i_j = nint(rj)

                    ii_s(ista) = i_i
                    jj_s(ista) = i_j

                    if(i_i .ge. 3 .and. i_i .le. ni-2 .and.
     1                 i_j .ge. 3 .and. i_j .le. nj-2            )then

                      if(trim(var_2d) .eq. 'SWI')then
                        if(iwrite .eq. iwrite/20*20)then
                          write(6,*)'sv '
                          write(6,*)'sv Sta   i    j   VIS frac tb8_k  '
     1                    //'t_gnd_k t_sfc_k cv_s_mx cvr_mx '
     1                    //'solalt 9pt  rad_fc '
     1                    //'rad_ob rad_th ratio cv_sol  df'
                        endif

!                       Calculate 9pt cover
                        cvr_9pt = 0.

!                       Calculate 25pt cover
                        cvr_25pt = 0.

                        iwrite = iwrite + 1

                        c1_c = ' '
                        c3_discrep = '   '

                        rad_ratio = 0.
                        cv_solar = 0.
                        cv_diff = 0.

                        rad2_s(ista) = var_s(ista)            
                        swi_s(ista) = swi_2d(i_i,i_j)
 
                        if(swi_s(ista) .gt. 0.)then
                          radob_ratio = var_s(ista) / swi_s(ista)         
                        else
                          radob_ratio = 1.0
                        endif
 
!                       QC checks
                        if(cvr_max(i_i,i_j) .le. .10 .and. 
     1                   radob_ratio .lt. 0.3      .and.
     1                   swi_s(ista) .ge. 100.           )then
                          c1_c = '-' ! Suspected low
                        endif

                        if(radob_ratio .gt. 3.0      .and.      
     1                     var_s(ista) .ge. 400.           )then
                           c1_c = '+' ! Suspected high
                        endif

                        if(radob_ratio .lt. 0.1 .and. 
     1                     swi_s(ista) .ge. 100.      )then
                          c1_c = '*' ! QC'd out
                          rad2_s(ista) = r_missing_data
                        endif

                        if(var_s(ista) - rad_clr(i_i,i_j) .gt. 500.)then       
                          if(rad_clr(i_i,i_j) .gt. 100.)then
                            if(var_s(ista)/rad_clr(i_i,i_j) .gt. 2.5
     1                                                             )then
                               c1_c = '*'
                               rad2_s(ista) = r_missing_data
                            endif
                          endif
                        endif

                        write(6,1111,err=1112)'sv',stations(ista)(1:3)
     1                             ,i_i,i_j
     1                             ,cloud_frac_vis_a(i_i,i_j)
     1                             ,tb8_k(i_i,i_j)
     1                             ,t_gnd_k(i_i,i_j)
     1                             ,t_sfc_k(i_i,i_j)
     1                             ,cvr_sao_max(i_i,i_j)
     1                             ,cvr_max(i_i,i_j)
     1                             ,solar_alt(i_i,i_j)
     1                             ,cvr_9pt
!    1                             ,cvr_25pt
     1                             ,swi_2d(i_i,i_j)
     1                             ,var_s(ista)
     1                             ,rad_clr(i_i,i_j)
     1                             ,rad_ratio
     1                             ,cv_solar
     1                             ,cv_diff
     1                             ,c1_c
1111                    format(1x,a2,1x,a3,2i5,f8.2,3f8.1,f7.2,f8.2,f6.1       
     1                        ,f6.2,f8.1,2f7.1,3f6.2,1x,a,1x)


                      elseif(trim(var_2d) .eq. 'TSF' .OR. 
     1                       trim(var_2d) .eq. 'DSF'     )then ! convert K to F
                        var_fcst_s(ista) = k_to_f(var_fcst_2d(i_i,i_j))       
                   
                      elseif(trim(var_2d) .eq. 'TPW')then      ! convert M to CM
                        var_fcst_s(ista) = var_fcst_2d(i_i,i_j) * 100.

                      elseif(trim(var_2d) .eq. 'R01' .OR. 
     1                       trim(var_2d) .eq. 'R03' .OR.
     1                       trim(var_2d) .eq. 'R06' .OR.
     1                       trim(var_2d) .eq. 'R24' .OR.
     1                       trim(var_2d) .eq. 'RTO'     )then ! convert M to IN

!                       if(var_fcst_2d(i_i,i_j) .gt. .000254 .OR. ! .01 inch
!    1                     var_s(ista)          .gt. 0.      )then
                        if(var_fcst_2d(i_i,i_j) .ge. 0.      .OR. ! .00 inch
     1                     var_s(ista)          .ge. 0.      )then
                            var_fcst_s(ista) = var_fcst_2d(i_i,i_j) 
     1                                       / .0254 
                           if(var_s(ista) .ge. 1.0)then
                              c1_c = '*'
                              write(6,1001)var_2d,ista
     1                                    ,stations(ista)(1:5),i_i,i_j        
     1                                    ,var_s(ista),var_fcst_s(ista)
     1                                    ,c1_c
                           endif
                        else
                            var_fcst_s(ista) = r_missing_data
                            var_s(ista) = r_missing_data
                        endif

                      else
                        var_fcst_s(ista) = var_fcst_2d(i_i,i_j)            

                      endif

                      if(var_fcst_s(ista)     .ne. r_missing_data .AND.
     1                   var_s(ista)          .ne. r_missing_data .AND.
     1                   var_fcst_2d(i_i,i_j) .ne. r_missing_data  )then   
                        if(cnt_write .le. 50.)then
                          c1_c = ' '
                          write(6,1001)var_2d,ista
     1                                ,stations(ista)(1:5),i_i,i_j 
     1                                ,var_s(ista),var_fcst_s(ista)
     1                                ,c1_c
1001                      format(1x,a3,' ob/fcst',i9,1x,a5,1x,2i8,2f9.3
     1                          ,1x,a1) 
                        endif
                        cnt_write = cnt_write + 1.
                      endif

                      sumobs = sumobs + var_s(ista)
                      sumanl = sumanl + swi_2d(i_i,i_j)
                      sumcld = sumcld + cvr_max(i_i,i_j)
                      sumalt = sumalt + solar_alt(i_i,i_j)
                      sumscl = sumscl + cvr_scl_a(i_i,i_j)
                      cnt = cnt + 1.

                      cvr_s(ista) = cvr_max(i_i,i_j)
1112                  continue

                    else ! outside domain
                      var_fcst_s(ista) = r_missing_data

                    endif ! ob is in domain

                  else
                    var_fcst_s(ista) = r_missing_data

                  endif ! valid ob                                  

                enddo ! ista

                write(6,*)
                write(6,*)' Generic stats, cnt = ',nint(cnt)

                if(cnt .eq. 0.)then
                  write(6,*)' No obs - skipping calculation of stats'
                  goto 980 ! Skip calculation of stats and write initialized flag values
                endif
 
1200            if(trim(var_2d) .eq. 'SWI')then
                  call stats_1d(maxsta,swi_s,rad2_s
     1                   ,'Solar Radiation (QCed): '
     1                   ,a_t,b_t,xbar,ybar
     1                   ,bias,std,r_missing_data,istatus)
                  rms_a(ifield,imodel,itime_fcst) = std
                  cnt_a(ifield,imodel,itime_fcst) = cnt
                elseif(trim(var_2d) .eq. 'TSF')then
                  call stats_1d(maxsta,var_fcst_s,var_s
     1                   ,'Surface Temperature     '
     1                   ,a_t,b_t,xbar,ybar
     1                   ,bias,std,r_missing_data,istatus)
                  rms_a(ifield,imodel,itime_fcst) = std
                  cnt_a(ifield,imodel,itime_fcst) = cnt
                elseif(trim(var_2d) .eq. 'DSF')then
                  call stats_1d(maxsta,var_fcst_s,var_s
     1                   ,'Surface Dewpoint     '
     1                   ,a_t,b_t,xbar,ybar
     1                   ,bias,std,r_missing_data,istatus)
                  rms_a(ifield,imodel,itime_fcst) = std
                  cnt_a(ifield,imodel,itime_fcst) = cnt
                elseif(trim(var_2d) .eq. 'USF')then
                  call stats_1d(maxsta,var_fcst_s,var_s
     1                   ,'Surface U Wind Component'
     1                   ,a_t,b_t,xbar,ybar
     1                   ,bias,std,r_missing_data,istatus)
                  rms_a(ifield,imodel,itime_fcst) = std
                  cnt_a(ifield,imodel,itime_fcst) = cnt
                elseif(trim(var_2d) .eq. 'VSF')then
                  call stats_1d(maxsta,var_fcst_s,var_s
     1                   ,'Surface V Wind Component'
     1                   ,a_t,b_t,xbar,ybar
     1                   ,bias,std,r_missing_data,istatus)
                  rms_a(ifield,imodel,itime_fcst) = std
                  cnt_a(ifield,imodel,itime_fcst) = cnt
                elseif(trim(var_2d) .eq. 'SSF')then
                  call stats_1d(maxsta,var_fcst_s,var_s
     1                   ,'Surface Wind Speed'
     1                   ,a_t,b_t,xbar,ybar
     1                   ,bias,std,r_missing_data,istatus)
                  rms_a(ifield,imodel,itime_fcst) = std
                  cnt_a(ifield,imodel,itime_fcst) = cnt
                elseif(trim(var_2d) .eq. 'WSF')then
                  threshval = -99.9
                  xbar = rms_a(4,imodel,itime_fcst) 
                  ybar = rms_a(5,imodel,itime_fcst)
                  if(xbar .ne. threshval .and. ybar .ne. threshval)then
                      std = sqrt(rms_a(4,imodel,itime_fcst)**2
     1                          +rms_a(5,imodel,itime_fcst)**2)
                  else
                      std = threshval
                  endif
                  rms_a(ifield,imodel,itime_fcst) = std
                  cnt_a(ifield,imodel,itime_fcst) = 
     1            cnt_a(4,imodel,itime_fcst)
                elseif(trim(var_2d) .eq. 'TPW')then
                  call stats_1d(maxsta,var_fcst_s,var_s
     1                   ,'Integrated Water Vapor'
     1                   ,a_t,b_t,xbar,ybar
     1                   ,bias,std,r_missing_data,istatus)
                  rms_a(ifield,imodel,itime_fcst) = std
                  cnt_a(ifield,imodel,itime_fcst) = cnt
                elseif(var_2d(1:2) .eq. 'R0' .OR. 
     1                 var_2d(1:2) .eq. 'R2'      )then      
                  call stats_1d(maxsta,var_fcst_s,var_s
     1                   ,'Incremental Precip'             
     1                   ,a_t,b_t,xbar,ybar
     1                   ,bias,std,r_missing_data,istatus)
                  rms_a(ifield,imodel,itime_fcst) = std
                  cnt_a(ifield,imodel,itime_fcst) = cnt
                elseif(trim(var_2d) .eq. 'RTO')then
                  call stats_1d(maxsta,var_fcst_s,var_s
     1                   ,'Run Total Precip'             
     1                   ,a_t,b_t,xbar,ybar
     1                   ,bias,std,r_missing_data,istatus)
                  rms_a(ifield,imodel,itime_fcst) = std
                  cnt_a(ifield,imodel,itime_fcst) = cnt
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
     1                           ,nint(cnt_a(ifield,imodel,itime_fcst))
!               write(39,710)atime_s,xbar,ybar,std
!               write(38,710)atime_s,xbar,ybar,std
710             format(1x,a24,3f10.3,i8)

              enddo ! itime_fcst

              write(lun_out,*)
              write(lun_out,*)

            endif ! c_model .ne. lga

          enddo ! model

          write(6,*)' Closing lun_out ',lun_out
          close (lun_out) 

        enddo ! fields

 999    write(6,*)' End of subroutine verif_fcst_pt_2d'

        return

        end
