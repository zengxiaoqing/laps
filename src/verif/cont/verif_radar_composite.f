          
        subroutine verif_radar_composite(i4time_sys,a9time,
     1                  model_fcst_intvl,
     1                  model_fcst_len,
     1                  model_cycle_time,
     1                  laps_cycle_time,
     1                  NX_L,NY_L,
     1                  NZ_L,
     1                  r_missing_data,
     1                  n_plot_times,
     1                  l_persist,
     1                  j_status)

        include 'lapsparms.for' ! maxbgmodels 

        real var_anal_3d(NX_L,NY_L,NZ_L)
        real var_fcst_3d(NX_L,NY_L,NZ_L)
        real rqc(NX_L,NY_L)

        logical lmask_and_3d(NX_L,NY_L,NZ_L)
        logical lmask_rqc_3d(NX_L,NY_L,NZ_L)
        logical l_col /.true./
        logical l_exist
        logical l_plot_criteria 
        logical l_persist 
        logical l_req_all_mdls / .true. /
        logical l_time_outcoord_hhmm 

!       integer       maxbgmodels
!       parameter     (maxbgmodels=10)
        character*30  c_fdda_mdl_src(maxbgmodels)
        character*30  c_fdda_mdl_hdr(maxbgmodels)
        character*256 cline
        character*1 char

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
        character*10 var_2d
        character*5 fcst_hh_mm
        character*9 a9time,a9time_valid,a9time_initial
        character*24 a24time_valid
        character*24 a24time_valid_file, a24time_valid_expected
        character*150 hist_dir, cont_dir, verif_dir, plot_dir
        character*150 hist_file, members_file
        character*150 bias_file_in, ets_file_in
        character*150 bias_file_out, ets_file_out
        character*150 summary_file_out
        character*10 compdir

        integer n_fields
        parameter (n_fields=6)
        character*10 var_a(n_fields)
        integer nthr_a(n_fields)   ! number of thresholds for each field
        integer istart_a(n_fields) ! start time for each field              
        character*2 c2_region
        character*10 c_thr

!       Specify what is being verified
        data var_a  /'REF','LMR','PCP_01','PCP_03','PCP_06','PCP_24'/ ! 3-D / composite ref
        data nthr_a   /5,5,7,7,7,7/        
        data istart_a /0,0,1,1,1,1/        

        integer,parameter :: k12 = selected_int_kind(12)
        integer (kind=k12) :: contable(0:1,0:1)
        integer (kind=k12) :: idenom_sum,idenom_sum2,ihits_sum                       

        integer maxthr
        parameter (maxthr=7)

        real pcp_thr(7)
        data pcp_thr /.01,.05,0.1,0.5,1.0,2.0,5.0/

        real cont_4d(NX_L,NY_L,NZ_L,maxthr)
        real bias(maxbgmodels,0:max_fcst_times,max_regions,maxthr)
        real ets(maxbgmodels,0:max_fcst_times,max_regions,maxthr)
        real bias_comp(maxbgmodels,0:max_fcst_times,max_regions,maxthr)
        real ets_comp(maxbgmodels,0:max_fcst_times,max_regions,maxthr)
        real 
     1  frac_coverage(maxbgmodels,0:max_fcst_times,max_regions,maxthr)
        real 
     1  frac_cvr_comp(maxbgmodels,0:max_fcst_times,max_regions,maxthr)
        integer 
     1  n(maxbgmodels,0:max_fcst_times,max_regions,maxthr,0:1,0:1)
        integer (kind=k12) :: 
     1  n_sum(maxbgmodels,0:max_fcst_times,max_regions,maxthr,0:1,0:1)

!       Diagnostic logging arrays
        real ets_max(maxbgmodels,0:max_fcst_times,maxthr)
        real ets_min(maxbgmodels,0:max_fcst_times,maxthr)
        real ets_sum(maxbgmodels,0:max_fcst_times,maxthr)
        real ets_count(maxbgmodels,0:max_fcst_times,maxthr)
        real ets_mean(maxbgmodels,0:max_fcst_times,maxthr)

!       integer nmissing_m(maxbgmodels)
        integer nsuccess_m(maxbgmodels)
        integer nincomplete_m(maxbgmodels)
        integer incomplete_run_m(maxbgmodels)
        integer n_plot_times_m(maxbgmodels,0:max_fcst_times,n_fields)

       write(6,*)
       write(6,*)' Start subroutine verif_radar_composite...'

       l_time_outcoord_hhmm = .true.

       do i_period = 1,2
        
        if(i_period .eq. 1)then
         ndays = 7
         compdir = 'comp'
        else
         ndays = 30
         compdir = 'comp2'
        endif

        write(6,*)' Processing stats for period/dir = ',ndays,compdir

        rmiss = -999.
        imiss = -999

        thresh_var = 20. ! lowest threshold for this variable

        write(6,*)' Time is ',i4time_sys,a9time
        write(6,*)' n_plot_times ',n_plot_times

        n_plot_times_m(:,:,:) = 1

        i4_initial = i4time_sys

!       Determine verification timing
        model_verif_intvl = max(laps_cycle_time,model_fcst_intvl)
        n_fcst_times = (model_fcst_len*60) / model_verif_intvl

        write(6,*)' model_verif_intvl = ',model_verif_intvl
        write(6,*)' n_fcst_times = ',n_fcst_times

        lun_in = 21

        lun_bias_in = 31
        lun_ets_in = 32
        lun_bias_out = 41
        lun_ets_out = 42
        lun_summary_out = 43

!       Get fdda_model_source and 'n_fdda_models' from static file
        call get_fdda_model_source(c_fdda_mdl_src,n_fdda_models,istatus)

        if(l_persist .eqv. .true.)then
            n_fdda_models = n_fdda_models + 1
            c_fdda_mdl_src(n_fdda_models) = 'persistence'
            write(6,*)' Adding persistence to fdda_models'
        endif

        write(6,*)' n_fdda_models = ',n_fdda_models
        write(6,*)' c_fdda_mdl_src = '
     1            ,(c_fdda_mdl_src(m),m=1,n_fdda_models)

!       Update array of n_plot_times_m for available model exceptions
        do imodel = 2,n_fdda_models

!           Advection available just for LMR out to 3 hours
            if(.false.)then
              if(trim(c_fdda_mdl_src(imodel)) .eq. 'advection')then
                do ifield = 1,n_fields 
                    if(var_a(ifield) .ne. 'LMR')then ! only LMR is available
                        n_plot_times_m(imodel,:,ifield) = 0
                    else ! block out forecasts > 3 hours
                        n_fcst_times_m = 10800 / model_verif_intvl
                        if(max_fcst_times .gt. n_fcst_times_m)then 
                            n_plot_times_m(imodel,
     1                                 n_fcst_times_m+1:max_fcst_times,
     1                                 ifield) = 0 
                        endif
                    endif
                enddo
              endif
            endif

!           Only LMR available from the HRRR, RR, RAP-NH
            if(trim(c_fdda_mdl_src(imodel)) .eq. 'wrf-hrrr' .OR.
     1         trim(c_fdda_mdl_src(imodel)) .eq. 'rr'       .OR.
     1         trim(c_fdda_mdl_src(imodel)) .eq. 'rap-nh'        )then
                do ifield = 1,n_fields ! only LMR is available
                    if(var_a(ifield) .ne. 'LMR')then
                        n_plot_times_m(imodel,:,ifield) = 0
                    endif
                enddo
            endif

!           Only LMR & PCP_03 available from the NAM, every 3 hours
            if(trim(c_fdda_mdl_src(imodel)) .eq. 'nam')then
                do ifield = 1,n_fields                                 
                    if(var_a(ifield) .eq. 'LMR'                          
     1            .OR. var_a(ifield) .eq. 'PCP_03')then ! LMR & PCP_03 avail                    
                        nam_fcst_intvl = 10800          ! every 3 hours
                        do itime_fcst = 0,n_fcst_times
                            i4_fcst = itime_fcst*model_verif_intvl      
                            if(i4_fcst .ne. 
     1                        (i4_fcst/nam_fcst_intvl)*nam_fcst_intvl
     1                                                             )then      
                                n_plot_times_m(imodel,itime_fcst,ifield)
     1                                                               = 0
                            endif
                        enddo
                    else 
                        n_plot_times_m(imodel,:,ifield) = 0
                    endif
                enddo
            endif

!           Only LMR & PCP_06 available from the NAM-NH, every 6 hours
            if(trim(c_fdda_mdl_src(imodel)) .eq. 'nam-nh')then
                do ifield = 1,n_fields                                 
                    if(var_a(ifield) .eq. 'LMR'                          
     1            .OR. var_a(ifield) .eq. 'PCP_06')then ! LMR & PCP_06 avail                    
                        nam_fcst_intvl = 21600          ! every 6 hours
                        do itime_fcst = 0,n_fcst_times
                            i4_fcst = itime_fcst*model_verif_intvl      
                            if(i4_fcst .ne. 
     1                        (i4_fcst/nam_fcst_intvl)*nam_fcst_intvl
     1                                                             )then      
                                n_plot_times_m(imodel,itime_fcst,ifield)
     1                                                               = 0
                            endif
                        enddo
                    else 
                        n_plot_times_m(imodel,:,ifield) = 0
                    endif
                enddo
            endif

!           Just use 1-hourly fields for PCP_01
            intvl_pcp = 3600
            do ifield = 1,n_fields                                 
                if(var_a(ifield) .eq. 'PCP_01')then ! PCP_01 avail every 1 hour
                    do itime_fcst = 0,n_fcst_times
                        i4_fcst = itime_fcst*model_verif_intvl      
                        if(i4_fcst .ne. (i4_fcst/intvl_pcp)*intvl_pcp
     1                                                             )then 
                            n_plot_times_m(imodel,itime_fcst,ifield)
     1                                                           = 0
                        endif
                    enddo
!               else 
!                   n_plot_times_m(imodel,:,ifield) = 0
                endif
            enddo

!           Just use 3-hourly fields for PCP_03
            intvl_pcp = 10800
            do ifield = 1,n_fields                                 
                if(var_a(ifield) .eq. 'PCP_03')then ! PCP_03 avail every 3 hours
                    do itime_fcst = 0,n_fcst_times
                        i4_fcst = itime_fcst*model_verif_intvl      
                        if(i4_fcst .ne. (i4_fcst/intvl_pcp)*intvl_pcp
     1                                                             )then 
                            n_plot_times_m(imodel,itime_fcst,ifield)
     1                                                           = 0
                        endif
                    enddo
!               else 
!                   n_plot_times_m(imodel,:,ifield) = 0
                endif
            enddo

!           Just use 6-hourly fields for PCP_06
            intvl_pcp = 21600
            do ifield = 1,n_fields                                 
                if(var_a(ifield) .eq. 'PCP_06')then ! PCP_06 avail every 6 hours
                    do itime_fcst = 0,n_fcst_times
                        i4_fcst = itime_fcst*model_verif_intvl      
                        if(i4_fcst .ne. (i4_fcst/intvl_pcp)*intvl_pcp
     1                                                             )then 
                            n_plot_times_m(imodel,itime_fcst,ifield)
     1                                                           = 0
                        endif
                    enddo
!               else 
!                   n_plot_times_m(imodel,:,ifield) = 0
                endif
            enddo

!           Just use 24-hourly fields for PCP_24
            intvl_pcp = 86400
            do ifield = 1,n_fields                                 
                if(var_a(ifield) .eq. 'PCP_24')then ! PCP_24 avail every 24 hrs
                    do itime_fcst = 0,n_fcst_times
                        i4_fcst = itime_fcst*model_verif_intvl      
                        if(i4_fcst .ne. (i4_fcst/intvl_pcp)*intvl_pcp
     1                                                             )then 
                            n_plot_times_m(imodel,itime_fcst,ifield)
     1                                                           = 0
                        endif
                    enddo
!               else 
!                   n_plot_times_m(imodel,:,ifield) = 0
                endif
            enddo

!           Only REF & LMR available from persistence
            if(trim(c_fdda_mdl_src(imodel)) .eq. 'persistence')then
                do ifield = 1,n_fields ! all unavailable except REF & LMR
                    if(var_a(ifield) .ne. 'REF' .AND.
     1                 var_a(ifield) .ne. 'LMR'       )then
                        n_plot_times_m(imodel,:,ifield) = 0
                    endif
                enddo
            endif

!           Assume WNDW-WARW isn't available for composite forecasts  
!           if(trim(c_fdda_mdl_src(imodel)) .eq. 'wndw-warw')then
!               n_plot_times_m(imodel,:,:) = 0
!           endif

!           Assume NAM isn't available for composite forecasts  
!           if(trim(c_fdda_mdl_src(imodel)) .eq. 'nam')then
!               n_plot_times_m(imodel,:,:) = 0
!           endif

!           Assume WRF-TOM isn't available for composite forecasts  
!           if(trim(c_fdda_mdl_src(imodel)) .eq. 'wrf-tom')then
!               n_plot_times_m(imodel,:,:) = 0
!           endif

        enddo ! imodel

!       Read in data file with region points
        n_models = n_fdda_models
!       call read_region_info(maxbgmodels,max_fcst_times,max_regions
!    1                       ,n_models,n_fcst_times,n_regions
!    1                       ,il,ih,jl,jh,lun_in)

!       In 'nest7grid.parms' model #1 is lga usually
!       In 'verif_regions.dat' model #1 represents the analysis

!       if(n_fdda_models .ne. n_models)then
!           write(6,*)' ERROR n_models differs from n_fdda_models '
!    1                       ,n_models,n_fdda_models
!           stop
!       endif

        call get_directory('verif',verif_dir,len_verif)

        n_init_times = ((ndays * 86400) / model_cycle_time) - 1

        do ifield = 1,n_fields

!        Initialize arrays
         bias = -999.
         ets = -999.
         frac_coverage = -999.
         n_sum = 0
        
!        Diagnostic array initialization
         ets_max   = -999.9
         ets_min   = +999.9
         ets_sum   =   0.
         ets_count =   0.

         var_2d = var_a(ifield)
         call s_len(var_2d,lenvar)

         write(6,*)' var / ndays / n_init_times = ',trim(var_2d)
     1                                             ,ndays,n_init_times       

         iregion = 1

         nmissing = 0
         nsuccess = 0
         nincomplete = 0

!        nmissing_m = 0
         nsuccess_m = 0
         nincomplete_m = 0

         frac_thr = 0.15
!        frac_thr = 0.035
         nmissing_thr = int((1. - frac_thr) * float(n_init_times+1))
         nsuccess_thr = (n_init_times+1) - nmissing_thr

         do init = 0,n_init_times

          incomplete_run = 0   ! based on any of the models
          incomplete_run_m = 0 ! based on each model

          n = imiss ! Initialize

          i4_initial = i4time_sys - (init * model_cycle_time)       
          call make_fnam_lp(i4_initial,a9time_initial,istatus)

          write(6,*)                                            
          write(6,*)' Processing model cycle ',a9time_initial

!         Read individual bias files
          write(6,*)' nthr before assigning = ',nthr
          nthr = nthr_a(ifield) ! nthr may be unset or have earlier been stepped on
          do idbz = 1,nthr

           if(var_2d(1:3) .ne. 'PCP')then
               rdbz = float(idbz*10) + 10
               write(c_thr,901)nint(rdbz)
 901           format(i2)
           else
               write(c_thr,902)nint(pcp_thr(idbz)*100.)
 902           format(i4.4)
           endif

           plot_dir = verif_dir(1:len_verif)//var_2d(1:lenvar)
     1                                      //'/plot'
!    1                                      //c_model(1:len_model)
!    1                                      //'/'
           len_plot = len_verif + 5 + lenvar ! + len_model

!          write GNUplot file for the time series of this model (region 1)
           bias_file_in = plot_dir(1:len_plot)//'/'
     1                                     //trim(c_thr)//'/'
     1                                     //a9time_initial     
     1                                     //'.bias'     

           write(6,*)'bias_file_in = ',bias_file_in

           ets_file_in  = plot_dir(1:len_plot)//'/'
     1                                     //trim(c_thr)//'/'
     1                                     //a9time_initial     
     1                                     //'.ets'      

           write(6,*)'ets_file_in = ',ets_file_in

           members_file  = verif_dir(1:len_verif)//'/'
     1                                     //'members.txt'      

           write(6,*)'members_file = ',members_file

           inquire(file=bias_file_in,exist=l_exist)
           if(.not. l_exist)then
               write(6,*)' WARNING: file does not exist:',bias_file_in
               goto958
           endif ! l_exist

           inquire(file=ets_file_in,exist=l_exist)
           if(.not. l_exist)then
               write(6,*)' ERROR: file does not exist:',ets_file_in       
               goto958
           endif ! l_exist

           open(lun_bias_in,file=bias_file_in,status='old')
           open(lun_ets_in,file=ets_file_in,status='old')

!          Read comment with model member names
           read(lun_bias_in,*,err=958,end=958)
           read(lun_ets_in,51,err=958,end=958) cline
 51        format(a)
           write(6,*)'cline = ',cline

           char = ' '
           call csplit(cline,c_fdda_mdl_hdr,nelems,n_fdda_models,char
     1                ,istatus)
           if(istatus .ne. 1)then
               write(6,*)' Read error in ets file header line...'
               goto958
           endif

           write(6,*)'c_fdda_mdl_hdr nelems/n_fdda_models = '
     1                              ,nelems,n_fdda_models

906        format('# ',30(1x,a))  

           if(l_col)then
!              Read bias and ets values
               do itime_fcst = 0,n_fcst_times
                   i4_valid = i4_initial + itime_fcst*model_verif_intvl      
                   call cv_i4tim_asc_lp(i4_valid,a24time_valid_expected
     1                                 ,istatus)
                   call left_justify(a24time_valid_expected)

                   read(lun_bias_in,911)a24time_valid,    
     1                 (bias(imodel,itime_fcst,iregion,idbz)
     1                              ,imodel=2,n_fdda_models)     
                   call left_justify(a24time_valid)
                   if(a24time_valid .ne.
     1                a24time_valid_expected)then
                       write(6,*)
     1                 ' WARNING: imodel / a24time (expected/file-1)'
     1                 ,imodel,itime,a24time_valid_expected            
     1                              ,a24time_valid              
                       goto958
                   endif

                   read(lun_ets_in,911)a24time_valid,    
     1                 (ets(imodel,itime_fcst,iregion,idbz)
     1                              ,imodel=2,n_fdda_models)     
                   call left_justify(a24time_valid)
                   if(a24time_valid .ne.
     1                a24time_valid_expected)then
                       write(6,*)
     1                 ' WARNING: imodel / a24time (expected/file-2)'
     1                 ,imodel,itime,a24time_valid_expected           
     1                              ,a24time_valid         
                       goto958
                   endif

911                format(a24,3x,20f12.3)
               enddo ! itime_fcst

!              Read N values from separate blocks                     
               do jn = 0,1
               do in = 0,1
                 read(lun_bias_in,*)
                 read(lun_bias_in,*)
                 do itime_fcst = 0,n_fcst_times
                   i4_valid = i4_initial + itime_fcst*model_verif_intvl      
                   call cv_i4tim_asc_lp(i4_valid,a24time_valid_expected
     1                                 ,istatus)
                   call left_justify(a24time_valid_expected)

                   read(lun_bias_in,51)cline
                   read(cline,912,err=913)a24time_valid,    
     1                 (n(imodel,itime_fcst,iregion,idbz,in,jn)
     1                              ,imodel=2,n_fdda_models)     
912                format(a24,3x,20i12.3)

                   if(n(imodel,itime_fcst,iregion,idbz,in,jn) 
     1                                            .lt. -1000)then
                       write(6,*)' WARNING: anomalous N value',
     1                     n(imodel,itime_fcst,iregion,idbz,in,jn)
                   endif

                   goto914
913                write(6,*)' read ERROR in bias file N vals '
                   write(6,*)' check if model config has changed'
                   write(6,*)in,jn,itime_fcst
                   write(6,*)cline
                   goto 955

914                continue

!                  Blank out N values if we're not expecting/desiring data
!                  where(n_plot_times_m(:,itime_fcst,ifield) 
!    1                                                 .eq. 0)
!                      n(:,itime_fcst,iregion,idbz,in,jn) = 0
!                  endwhere

!                  write(6,915)in,jn,itime_fcst,
!    1                 (n(imodel,itime_fcst,iregion,idbz,in,jn)
!    1                              ,imodel=2,n_fdda_models)     
915                format('in,jn,itime,n',i2,i2,i4,20i10)

                   call left_justify(a24time_valid)
                   if(a24time_valid .ne.
     1                a24time_valid_expected)then
                       write(6,*)
     1                 ' WARNING: imodel / a24time (expected/file-3)'
     1                 ,imodel,itime,a24time_valid_expected           
     1                              ,a24time_valid              
                       goto958
                   endif
                 enddo ! itime_fcst
               enddo ! jn
               enddo ! in

!              Test for missing data in all times/models for this threshold
               do itime_fcst = istart_a(ifield),n_fcst_times
                 do imodel = 2,n_fdda_models
                   i_good_timestep_model = 0
                   do in = 0,1
                   do jn = 0,1
                     if(n(imodel,itime_fcst,iregion,idbz,in,jn)
     1                                                  .gt. 0)then
                       i_good_timestep_model = 1
                     endif
                   enddo ! jn
                   enddo ! in

!                  Flag as incomplete if missing during expected time
                   if(i_good_timestep_model .eq. 0 .AND.
     1                n_plot_times_m(imodel,itime_fcst,ifield) .eq. 1 
     1                                                   )then

                     if(incomplete_run .eq. 0)then       
                       write(6,916)init,itime_fcst,a9time_initial,imodel
916                    format(
     1                 ' WARNING: missing N values for init/time/model '
     1                        ,2i6,1x,a9,1x,i6)                      
                       incomplete_run = 1 
                     endif

                     if(incomplete_run_m(imodel) .eq. 0)then       
                       incomplete_run_m(imodel) = 1
                     endif

                   endif ! i_good_timestep_model = 0

                 enddo ! im
               enddo ! itime_fcst

!              Read fractional coverage values
               read(lun_bias_in,*)
               read(lun_bias_in,*)
               do itime_fcst = 0,n_fcst_times
                   i4_valid = i4_initial + itime_fcst*model_verif_intvl      
                   call cv_i4tim_asc_lp(i4_valid,a24time_valid_expected
     1                                 ,istatus)
                   call left_justify(a24time_valid_expected)

                   read(lun_bias_in,923,err=925)a24time_valid,    
     1                 (frac_coverage                
     1                   (imodel,itime_fcst,iregion,idbz)
     1                              ,imodel=2,n_fdda_models)     
923                format(a24,3x,20f12.5)
925                continue               
                   call left_justify(a24time_valid)
                   if(a24time_valid .ne.
     1                a24time_valid_expected)then
                       write(6,*)
     1                 ' WARNING: imodel / a24time (expected/file-4)'
     1                 ,imodel,itime,a24time_valid_expected           
     1                              ,a24time_valid                           
                       goto958
                   endif
               enddo ! itime_fcst

!              Read members.txt file
               open(lun_mem,file=members_file,status='old')
               do imodel = 2,n_fdda_models
                 read(lun_mem,*)c_model ! from members.txt               

                 write(6,*)' c_model check (members/namelist/header): '
     1                                       ,imodel,c_model
     1                                       ,c_fdda_mdl_src(imodel)
     1                                       ,c_fdda_mdl_hdr(imodel)

                 if(trim(c_model) .ne. 'persistence')then

!                  Compare members.txt file and namelist fdda parms
                   if(trim(c_model) .ne. 
     1                              trim(c_fdda_mdl_src(imodel)))then
                       write(6,*)
     1             ' ERROR: Models did not match (members.txt/namelist)'
                       close(lun_mem)
                       return
                   endif

!                  Compare ets file header and namelist fdda parms
                   if(trim(c_model) .ne. 
     1                              trim(c_fdda_mdl_hdr(imodel)))then
                       write(6,*)
     1       ' WARNING: Models did not match (ets file header/namelist)'
                       close(lun_mem)
                       goto 958
                   endif

                 endif ! other than persistence from members.txt

               enddo ! imodel
               close(lun_mem)

           endif

940        close(lun_bias_in)
           close(lun_ets_in)

          enddo ! idbz 

!         Update arrays for all models
!         where(n_plot_times_m(:,ifield) .gt. -1)
              nsuccess_m(:) = nsuccess_m(:) + (1 - incomplete_run_m(:))
!         endwhere

          nincomplete = nincomplete + incomplete_run
          nincomplete_m(:) = nincomplete_m(:) + incomplete_run_m(:)

          nincomplete_t = 0
          do imodel=2,n_fdda_models
              nincomplete_t = nincomplete_t + incomplete_run_m(imodel)
          enddo 

          if(nincomplete_t .eq. 0 .OR. 
     1                             (l_req_all_mdls .eqv. .false.) )then
              write(6,*)' Accumulating sums for this run '
              nsuccess = nsuccess + 1

              where (n(:,:,:,:,:,:) .ne. imiss)                
               n_sum(:,:,:,:,:,:) = n_sum(:,:,:,:,:,:) + n(:,:,:,:,:,:)
              end where

!             Write to log for diagnostic purposes
              if(var_2d .eq. 'LMR')then
                 write(6,*)'...............'

                 itime_start =  3600/model_verif_intvl
                 itime_end   = 21600/model_verif_intvl
                 do itime_fcst = itime_start,itime_end 
     1                          ,itime_end-itime_start 

                   write(6,*)'...............'

                   do imodel=2,n_fdda_models

!                    Calculated summed skill scores so far
                     contable(:,:) = 
     1                   n_sum(imodel,itime_fcst,1,1,:,:)

                     idbz = 1

                     call skill_scores(contable,0                      ! I
     1                                ,frac_cvr_val                    ! O
     1                                ,frac_obs                        ! O
     1                                ,frac_fcst                       ! O
     1                                ,bias_val                        ! O
     1                                ,ets_val)                        ! O

                     ihits_sum  = contable(0,0)
                     idenom_sum = contable(1,0) + contable(0,1) 
     1                          + contable(0,0) 

!                    idenom_sum2 = frac_cvr_val * float(NX_L*NY_L) 
!    1                           * ets_count(imodel,idbz)

                     ets_min(imodel,itime_fcst,idbz) = 
     1                   min(ets_min(imodel,itime_fcst,idbz),ets_val)

                     ets_max(imodel,itime_fcst,idbz) = 
     1                   max(ets_max(imodel,itime_fcst,idbz),ets_val)

                     if(ets_val .ne. rmiss .and. 
     1                  ets_min(imodel,itime_fcst,idbz) .ne. rmiss) then ! valid ets 
                         ets_sum(imodel,itime_fcst,idbz) = 
     1                   ets_sum(imodel,itime_fcst,idbz) + ets_val
                         ets_count(imodel,itime_fcst,idbz) = 
     1                   ets_count(imodel,itime_fcst,idbz) + 1.0
                         ets_mean(imodel,itime_fcst,idbz) =  
     1                   ets_sum(imodel,itime_fcst,idbz) / 
     1                   ets_count(imodel,itime_fcst,idbz)       
                     else
                         ets_mean(imodel,itime_fcst,idbz) = rmiss
                     endif

                     idenom_sum2 = frac_cvr_val * float(NX_L*NY_L) 
     1                           * ets_count(imodel,itime_fcst,idbz)

!                    Process contingency table for run
                     contable(:,:) = 
     1                   n(imodel,itime_fcst,1,1,:,:)

                     idenom_run = contable(1,0) + contable(0,1) 
     1                          + contable(0,0) 

!                    if(frac_coverage(imodel,itime_fcst,iregion,idbz)
!    1                  .ne. -999.)then
!                        idenom_run2 = 
!    1                   frac_coverage(imodel,itime_fcst,iregion,idbz)
!    1                    * float(NX_L*NY_L)
!                    else
!                        idenom_run2 = -999
!                    endif

                     write(6,950)a9time_initial,c_fdda_mdl_src(imodel)
     1                    ,var_2d
     1                    ,itime_fcst
     1                    ,n(imodel,itime_fcst,1,1,1,1)     ! negs run
     1                    ,n_sum(imodel,itime_fcst,1,1,1,1) ! negs sum     
     1                    ,n(imodel,itime_fcst,1,1,0,0)     ! hits run
     1                    ,n_sum(imodel,itime_fcst,1,1,0,0) ! hits sum     
     1                    ,idenom_run,idenom_sum            ! denoms               
     1                    ,nint(ets_count(imodel,itime_fcst,idbz))
     1                    ,ets(imodel,itime_fcst,iregion,idbz)       
     1                    ,ets_val
     1                    ,ets_mean(imodel,itime_fcst,idbz)
     1                    ,ets_min(imodel,itime_fcst,idbz)
     1                    ,ets_max(imodel,itime_fcst,idbz)

950                  format(a9,1x,a10,1x,a3,' itime_fcst',i3
     1                      ,' Negs (run/sum):',2i11
     1                      ,' Hits (run/sum):',2i10
     1                      ,' Denom (run/sum):',2i10 
     1                      ,' ETS (run/sum/mean/min/max) = '
     1                      ,i3,3f9.3,2x,2f9.3)

                   enddo ! imodel
                 enddo ! itime_fcst
              endif ! var_2d = 'LMR'
          else
              write(6,*)' Not accumulating sums for this run'

          endif ! accumulate sums for this run

          write(6,953)a9time_initial,nincomplete_t,
     1                (incomplete_run_m(imodel),imodel=2,n_fdda_models)
953       format(' incomplete_run_m at ',a9,' is ',i3,4x,20i3)   

955       close(lun_bias_in)                                      
          close(lun_ets_in)                                         
      
          goto 960 ! success for this time

!         Error Condition for this time
958       nmissing = nmissing + 1

960      enddo                    ! init (initialization time)

         write(6,*)'init neg ',n_sum(2,0,1,1,1,1)

         write(6,965)nsuccess,n_init_times+1,nsuccess_thr
     1          ,((float(nsuccess) / float(n_init_times+1))) * 100.
965      format(' success count ',i4,' out of ',i4,i4,' were needed' 
     1         ,f8.2,'%')
         write(6,966)(nsuccess_m(imodel),imodel=2,n_fdda_models)
966      format(' nsuccess_m is ',20i5)
         write(6,967)                        
     1      ( ((float(nsuccess_m(imodel)) / float(n_init_times+1))*100.)
     1                                        ,imodel=2,n_fdda_models) 
967      format(' nsuccess_m % is ',20f7.2)

         write(6,*)'nmissing is ',nmissing

         write(6,*)'nincomplete is ',nincomplete
         write(6,968)(nincomplete_m(imodel),imodel=2,n_fdda_models)
968      format(' nincomplete_m is ',20i5)   

         nruns_plotted = 0
         do imodel=2,n_fdda_models
             if(nsuccess_m(imodel) .ge. nsuccess_thr)then 
                 nruns_plotted = nruns_plotted + 1
             endif
         enddo 

         if(nsuccess .lt. nsuccess_thr)then
             l_plot_criteria = .false.
         else
             l_plot_criteria = .true.
         endif

         write(6,*)'l_plot_criteria = ',l_plot_criteria

!        Calculate composite bias/ets
         do idbz = 1,nthr
           do imodel=2,n_fdda_models
             do itime_fcst = 0,n_fcst_times
                 contable(0,0) = 
     1                 n_sum(imodel,itime_fcst,iregion,idbz,0,0)
                 contable(1,0) = 
     1                 n_sum(imodel,itime_fcst,iregion,idbz,1,0)
                 contable(0,1) = 
     1                 n_sum(imodel,itime_fcst,iregion,idbz,0,1)
                 contable(1,1) = 
     1                 n_sum(imodel,itime_fcst,iregion,idbz,1,1)

                 write(6,*)'init neg ',n_sum(2,0,1,1,1,1)
     1                ,n_sum(imodel,itime_fcst,iregion,idbz,1,1)
     1                ,imodel,itime_fcst,iregion,idbz,1,1           

                 write(6,970)c_fdda_mdl_src(imodel),itime_fcst
     1                      ,contable
!    1                      ,contable(0,0),contable(1,0)   
!    1                      ,contable(0,1),contable(1,1)   
970              format(/' Calling skill scores for ',a10,i3,3x,8i10)       

!                Test whether this model satisfies completeness criteria
!                Plots will show up for each model that has thresh % of runs with a complete set of forecast times
                 if(nsuccess_m(imodel) .ge. nsuccess_thr)then ! satisfies completeness criteria
                   lun_out = 6
                   call skill_scores(contable,lun_out                  ! I
     1                  ,frac_cvr_comp(imodel,itime_fcst,iregion,idbz) ! O
     1                  ,frac_obs_comp                                 ! O
     1                  ,frac_fcst                                     ! O
     1                  ,bias_comp(imodel,itime_fcst,iregion,idbz)     ! O
     1                  , ets_comp(imodel,itime_fcst,iregion,idbz))    ! O

                 else                                         ! does not satisfy criteria
                   frac_cvr_comp(imodel,itime_fcst,iregion,idbz) = rmiss
                   bias_comp(imodel,itime_fcst,iregion,idbz) = rmiss
                   ets_comp(imodel,itime_fcst,iregion,idbz) = rmiss

                 endif

             enddo ! itime_fcst
           enddo ! imodel
         enddo ! idbz

!        Define and write to summary*.txt file
         summary_file_out = verif_dir(1:len_verif)//var_2d(1:lenvar)
     1                              //'/plot'
     1                              //'/summary_'//trim(compdir)//'.txt'       

         write(6,*)'summary_file_out = ',summary_file_out

         open(lun_summary_out,file=summary_file_out,status='unknown')
         do imodel=2,n_fdda_models
             ipct = nint(  (float(nsuccess_m(imodel)) 
     1                    / float(n_init_times+1))*100.)
             write(lun_summary_out,969)ipct                        
 969         format(i3)
         enddo ! imodel
         write(lun_summary_out,*)l_plot_criteria
         ipct = nint(  (float(nsuccess          ) 
     1                / float(n_init_times+1))*100.)
         write(lun_summary_out,969)ipct                        
!        write(lun_summary_out,*)frac_obs_comp
         close(lun_summary_out)

         if(nsuccess .lt. nsuccess_thr)then
             write(6,*)' Insufficient successful times to plot'
             goto 980
         endif

         write(6,*)
         write(6,*)
     1        ' Write output composite bias/ets files ############# '                                           
         write(6,*)' nthr before assigning = ',nthr
         nthr = nthr_a(ifield) ! nthr may be unset or have earlier been stepped on
         do idbz = 1,nthr

           if(var_2d(1:3) .ne. 'PCP')then
               rdbz = float(idbz*10) + 10
               write(c_thr,901)nint(rdbz)
           else
               write(c_thr,902)nint(pcp_thr(idbz)*100.)
           endif
 
           plot_dir = verif_dir(1:len_verif)//var_2d(1:lenvar)
     1                                      //'/plot'
!    1                                      //c_model(1:len_model)
!    1                                      //'/'
           len_plot = len_verif + 5 + lenvar ! + len_model

           i4_initial = i4time_sys                                   
           call make_fnam_lp(i4_initial,a9time_initial,istatus)

!          write GNUplot file for the time series of this model (region 1)
           bias_file_out = plot_dir(1:len_plot)//'/'
     1                                     //trim(c_thr)
     1                                     //'_'//trim(compdir)//'/'
     1                                     //a9time_initial     
     1                                     //'.bias'     

           write(6,*)'bias_file_out = ',bias_file_out

           ets_file_out  = plot_dir(1:len_plot)//'/'
     1                                     //trim(c_thr)
     1                                     //'_'//trim(compdir)//'/'
     1                                     //a9time_initial     
     1                                     //'.ets'      

           write(6,*)'ets_file_out = ',ets_file_out

           open(lun_bias_out,file=bias_file_out,status='unknown')
           open(lun_ets_out,file=ets_file_out,status='unknown')

!          Write comment with model member names
           write(lun_bias_out,906)(trim(c_fdda_mdl_src(imodel))
     1                              ,imodel=2,n_fdda_models)  
           write(lun_ets_out,906)(trim(c_fdda_mdl_src(imodel))
     1                              ,imodel=2,n_fdda_models)  
           if(l_col)then
!              Write bias and ets values
               do itime_fcst = 0,n_fcst_times
                   i4_fcst = itime_fcst*model_verif_intvl      
                   if(l_time_outcoord_hhmm .eqv. .true.)then
                       i4_fcst_hh = i4_fcst / 3600
                       i4_fcst_mm = i4_fcst/60 - (i4_fcst_hh*60)
                       write(a24time_valid,975)i4_fcst_hh,i4_fcst_mm
 975                   format(i2.2,1x,i2.2,19x)
                   else
                       i4_valid = i4_initial + i4_fcst
                       call cv_i4tim_asc_lp(i4_valid,a24time_valid
     1                                     ,istatus)
                   endif

                   write(lun_bias_out,911)a24time_valid,    
     1                 (bias_comp(imodel,itime_fcst,iregion,idbz)
     1                              ,imodel=2,n_fdda_models)     
                   write(lun_ets_out,911)a24time_valid,    
     1                 (ets_comp(imodel,itime_fcst,iregion,idbz)
     1                              ,imodel=2,n_fdda_models)     
               enddo ! itime_fcst

!              Write n values in separate blocks                     
               do jn = 0,1
               do in = 0,1
                 write(lun_bias_out,*)
                 write(lun_bias_out,*)
                 do itime_fcst = 0,n_fcst_times
                   i4_fcst = itime_fcst*model_verif_intvl      
                   if(l_time_outcoord_hhmm .eqv. .true.)then
                       i4_fcst_hh = i4_fcst / 3600
                       i4_fcst_mm = i4_fcst/60 - (i4_fcst_hh*60)
                       write(a24time_valid,975)i4_fcst_hh,i4_fcst_mm
                   else
                       i4_valid = i4_initial + i4_fcst
                       call cv_i4tim_asc_lp(i4_valid,a24time_valid
     1                                     ,istatus)
                   endif

                   write(6           ,912)a24time_valid,    
     1                 (n_sum(imodel,itime_fcst,iregion,idbz,in,jn)
     1                              ,imodel=2,n_fdda_models)     
                   write(lun_bias_out,912)a24time_valid,    
     1                 (n_sum(imodel,itime_fcst,iregion,idbz,in,jn)
     1                              ,imodel=2,n_fdda_models)     
                 enddo ! itime_fcst
               enddo ! jn
               enddo ! in

!              Write fractional coverage values
               write(lun_bias_out,*)
               write(lun_bias_out,*)
               do itime_fcst = 0,n_fcst_times
                   i4_fcst = itime_fcst*model_verif_intvl      
                   if(l_time_outcoord_hhmm .eqv. .true.)then
                       i4_fcst_hh = i4_fcst / 3600
                       i4_fcst_mm = i4_fcst/60 - (i4_fcst_hh*60)
                       write(a24time_valid,975)i4_fcst_hh,i4_fcst_mm
                   else
                       i4_valid = i4_initial + i4_fcst
                       call cv_i4tim_asc_lp(i4_valid,a24time_valid
     1                                     ,istatus)
                   endif

                   write(lun_bias_out,923)a24time_valid,    
     1                 (frac_cvr_comp                
     1                   (imodel,itime_fcst,iregion,idbz)
     1                              ,imodel=2,n_fdda_models)     
               enddo ! itime_fcst

           endif

           close(lun_bias_out)
           close(lun_ets_out)

         enddo ! idbz

 980    enddo ! ifield 

       enddo ! i_period

 999   write(6,*)' End of subroutine verif_radar_composite'

       return

       end

