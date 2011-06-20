
        program verif_fcst_main

        character*9 a9time
	character*300 dir_t,filenamet
        character*10 c_n_fcst_times,c_model_cycle_time_sec
        integer, parameter :: lun=120

!       This can be changed to read the modeltime for the forecast
!       initialization?
!       call get_systime(i4time,a9time,istatus)
!       if(istatus .ne. 1)go to 999

        ISTAT = init_timer()

        call get_directory('time',dir_t,istatus)
        call s_len(dir_t,len_dir_t)
        filenamet = dir_t(1:len_dir_t)//'/modeltime.dat'
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

        maxsta = 1000000

!       Use getarg for model cycle time
        call getarg(1,c_model_cycle_time_sec)
        read(c_model_cycle_time_sec,*)model_cycle_time_sec

!       model_cycle_time_sec = 3600  
!       model_cycle_time_sec = 10800

!       Use getarg for number of forecasts
        call getarg(2,c_n_fcst_times)
        read(c_n_fcst_times,*)n_fcst_times

        write(6,*)' model_cycle_time_sec/n_fcst_times (from getarg) = '  
     1             ,model_cycle_time_sec,n_fcst_times
          
        call verif_fcst_pt_2d(i4time,a9time,laps_cycle_time,
     1                     NX_L,NY_L,
     1                     NZ_L,
     1                     maxsta,
     1                     r_missing_data,
     1                     model_cycle_time_sec,
     1                     n_fcst_times,
     1                     j_status)

        call verif_fcst_pt_3d(i4time,a9time,laps_cycle_time,
     1                     NX_L,NY_L,
     1                     NZ_L,
     1                     maxsta,
     1                     r_missing_data,
     1                     model_cycle_time_sec,
     1                     n_fcst_times,
     1                     j_status)

999     continue

        end
          
        subroutine verif_fcst_pt_2d(i4time_sys,a9time,laps_cycle_time,
     1                  ni,nj,
     1                  nk,
     1                  maxsta,
     1                  r_missing_data,
     1                  model_cycle_time_sec,
     1                  n_fcst_times,
     1                  j_status)

        include 'read_sfc.inc'

        real var_anal_2d(ni,nj)
        real var_fcst_2d(ni,nj)
        real lat(ni,nj)
        real lon(ni,nj)
        real topo(ni,nj)
        real rlaps_land_frac(ni,nj)

        real k_to_f

        logical lmask_and_3d(ni,nj,nk)
        logical lmask_rqc_3d(ni,nj,nk)

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
        character*3 var_2d
        character*9 a9time,a9time_valid,a9time_init
        character*24 atime_s
        character*150 hist_dir, cont_dir, verif_dir
        character*150 hist_file

        integer n_fields
        parameter (n_fields=2)
        character*10 ext_anal_a(n_fields), ext_fcst_a(n_fields)
        character*10 var_a(n_fields)
        character*2 c2_region

!       Specify what is being verified
        data ext_fcst_a /'fsf','fsf'/ ! 2-D
!       data ext_anal_a /'lps'/ ! 3-D reflectivity
        data var_a      /'SWI','TSF'/ ! 2-D downward solar radiation

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

        real cld_snd(max_cld_snd,KCLOUD)
        integer ista_snd(max_cld_snd)
        real cld_hts(KCLOUD)

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

        write(6,*)' Start subroutine verif_fcst_pt_2d...'

        call get_laps_domain_95(ni,nj,lat,lon,topo
     1           ,rlaps_land_frac,grid_spacing_cen_m,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error getting LAPS domain'
            return
        endif

!       n_fcst_times = 2 ! 38

        cnt = 0.

        thresh_var = 20. ! lowest threshold for this variable

        i4_initial = i4time_sys

        lun_in = 21

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

          if(trim(var_2d) .eq. 'SWI')then
              istart = 1
          else
              istart = 0
          endif

          lun_out = 39

          hist_dir = verif_dir(1:len_verif)//var_2d(1:lenvar)
     1                                       //'/pt'
!    1                                       //c_model(1:len_model)
          len_hist = len_verif + 3 + lenvar ! + len_model

          cont_dir = verif_dir(1:len_verif)//var_2d(1:lenvar)
     1                                       //'/cont/'
     1                                       //c_model(1:len_model)
     1                                       //'/'
          len_cont = len_verif + 6 + lenvar + len_model

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

                write(6,*)' Reading surface obs - i4_valid = ',i4_valid

!               Read surface obs
                call read_surface_data(i4_valid,atime_s,n_obs_g,n_obs_b, !regular LSO
     &         obstime,wmoid,stations,provider,wx_s,reptype,autostntype,
     &         lat_s,lon_s,elev_s,t_s,td_s,rh_s,dd_s,ff_s,ddg_s,ffg_s,
     &         alt_s,pstn_s,pmsl_s,delpch,delp,vis_s,solar_s,sfct,sfcm,
     &         pcp1,pcp3,pcp6,pcp24,snow,kloud_s,max24t,min24t,t_ea,
     &         td_ea,rh_ea,dd_ea,ff_ea,alt_ea,p_ea,vis_ea,solar_ea,
     &         sfct_ea,sfcm_ea,pcp_ea,snow_ea,store_amt,store_hgt,
     &         maxsta,jstatus)

                write(6,*)' Number of obs in box  ',n_obs_b,atime_s

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

                endif ! .true.

                if(var_2d .eq. 'SWI')then
                  swi_2d = var_fcst_2d
                endif

                do ista = 1,maxsta  
!                 stn(ista) = c_stations(ista)(1:3)
                  stn(ista) = stations(ista)(1:3)

!                 write(6,*)ista,var_s(ista)
                  swi_s(ista) = r_missing_data 
                  var_fcst_s(ista) = r_missing_data 
                  rad2_s(ista) = r_missing_data 

                  if(var_s(ista) .gt. threshval)then 
                    call latlon_to_rlapsgrid(lat_s(ista),lon_s(ista),lat
     1                          ,lon,ni,nj,ri,rj,istatus)

                    i_i = nint(ri)
                    i_j = nint(rj)

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

!                     Calculate 9pt cover
                      cvr_9pt = 0.
                      do ic = -1,1
                      do jc = -1,1
                        cvr_9pt = cvr_9pt + cvr_max(i_i+ic,i_j+jc)
                      enddo ! jc
                      enddo ! ic
                      cvr_9pt = cvr_9pt / 9.

!                     Calculate 25pt cover
                      cvr_25pt = 0.
                      do ic = -2,2
                      do jc = -2,2
                        cvr_25pt = cvr_25pt + cvr_max(i_i+ic,i_j+jc)
                      enddo ! jc
                      enddo ! ic
                      cvr_25pt = cvr_25pt / 25.

                      iwrite = iwrite + 1

                      if(var_2d .eq. 'SWI')then
                        c1_c = ' '
                        c3_discrep = '   '

                        if(solar_alt(i_i,i_j) .gt. 0.)then
                          rad_ratio = var_s(ista) / rad_clr(i_i,i_j)
                          cv_solar = (1.0 - rad_ratio)       ! 100% cloud cover 
     1                             / (cvr_scl_a(i_i,i_j))    ! has (1-cvr_scl) of possible
                                                             ! solar radiation
                          cv_diff = cv_solar - cvr_max(i_i,i_j)

                        else
                          rad_ratio = 0.
                          cv_solar = 0.
                          cv_diff = 0.

                        endif

                        rad2_s(ista) = var_s(ista)            
                        swi_s(ista) = swi_2d(i_i,i_j)
 
                        if(swi_s(ista) .gt. 0.)then
                          radob_ratio = var_s(ista) / swi_s(ista)         
                        else
                          radob_ratio = 1.0
                        endif

                      elseif(var_2d .eq. 'TSF')then
                        var_fcst_s(ista) = k_to_f(var_fcst_2d(i_i,i_j))       

                      endif
 
!                     QC checks
                      if(cvr_max(i_i,i_j) .le. .10 .and. 
     1                 radob_ratio .lt. 0.3      .and.
     1                 swi_s(ista) .ge. 100.           )then
                        c1_c = '-' ! Suspected low
                      endif

                      if(radob_ratio .gt. 3.0      .and.      
     1                   var_s(ista) .ge. 400.           )then
                         c1_c = '+' ! Suspected high
                      endif

                      if(radob_ratio .lt. 0.1 .and. 
     1                   swi_s(ista) .ge. 100.      )then
                        c1_c = '*' ! QC'd out
                        rad2_s(ista) = r_missing_data
                      endif

                      if(var_s(ista) - rad_clr(i_i,i_j) .gt. 500.)then
                        if(rad_clr(i_i,i_j) .gt. 100.)then
                          if(var_s(ista)/rad_clr(i_i,i_j) .gt. 2.5)then
                             c1_c = '*'
                             rad2_s(ista) = r_missing_data
                          endif
                        endif
                      endif

                      write(6,1111,err=1112)'sv',stations(ista)(1:3)
     1                           ,i_i,i_j
     1                           ,cloud_frac_vis_a(i_i,i_j)
     1                           ,tb8_k(i_i,i_j)
     1                           ,t_gnd_k(i_i,i_j)
     1                           ,t_sfc_k(i_i,i_j)
     1                           ,cvr_sao_max(i_i,i_j)
     1                           ,cvr_max(i_i,i_j)
     1                           ,solar_alt(i_i,i_j)
     1                           ,cvr_9pt
!    1                           ,cvr_25pt
     1                           ,swi_2d(i_i,i_j)
     1                           ,var_s(ista)
     1                           ,rad_clr(i_i,i_j)
     1                           ,rad_ratio
     1                           ,cv_solar
     1                           ,cv_diff
     1                           ,c1_c
1111                  format(1x,a2,1x,a3,2i5,f8.2,3f8.1,f7.2,f8.2,f6.1       
     1                      ,f6.2,f8.1,2f7.1,3f6.2,1x,a,1x)

                      sumobs = sumobs + var_s(ista)
                      sumanl = sumanl + swi_2d(i_i,i_j)
                      sumcld = sumcld + cvr_max(i_i,i_j)
                      sumalt = sumalt + solar_alt(i_i,i_j)
                      sumscl = sumscl + cvr_scl_a(i_i,i_j)
                      cnt = cnt + 1.

                      cvr_s(ista) = cvr_max(i_i,i_j)

1112                endif ! ob is in domain
                  endif ! ista .ne. 0 (valid value)
                enddo ! ista

                write(6,*)
                write(6,*)' Generic stats, cnt = ',nint(cnt)
                if(var_2d .eq. 'SWI')then
                  call stats_1d(maxsta,swi_s,rad2_s
     1                   ,'Solar Radiation (QCed): '
     1                   ,a_t,b_t,xbar,ybar
     1                   ,bias,std,r_missing_data,istatus)
                elseif(var_2d .eq. 'TSF')then
                  call stats_1d(maxsta,var_fcst_s,var_s
     1                   ,'Surface Temperature     '
     1                   ,a_t,b_t,xbar,ybar
     1                   ,bias,std,r_missing_data,istatus)
                endif

                write(6,*)
                write(6,*)' Writing to lun_out ',lun_out
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

 999    write(6,*)' End of subroutine verif_fcst_pt_2d'

        return

        end
