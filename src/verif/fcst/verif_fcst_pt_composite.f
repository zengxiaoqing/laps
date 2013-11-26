! Hongli Jiang: NX_L, NY_L, and NZ_L were commented out. Uncomment them
! and also added in verif_fcst where verif_fcst_compostie is called. 10/14/2013          

        subroutine verif_fcst_composite(i4time_sys,a9time,
     1                  model_fcst_intvl,
     1                  model_fcst_len,
     1                  model_cycle_time_in,
     1                  laps_cycle_time,
     1                  NX_L,NY_L,
     1                  NZ_L,
     1                  r_missing_data,
     1                  n_plot_times,
     1                  j_status)

        include 'lapsparms.for' ! maxbgmodels

        logical l_col /.true./
        logical l_exist
        logical l_plot_criteria
        logical l_req_all_mdls /.true./
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
        character*3 var_2d
        character*5 fcst_hh_mm
        character*9 a9time,a9time_valid,a9time_initial
        character*24 a24time_valid                                                
        character*24 a24time_valid_file, a24time_valid_expected
        character*150 hist_dir, cont_dir, verif_dir, plot_dir
        character*150 hist_file, members_file
        character*150 stats_file_in                 
        character*150 stats_file_out                 
        character*150 summary_file_out
        character*10 compdir

        integer n_fields
        parameter (n_fields=13)
        character*10 var_a(n_fields)
        integer nthr_a(n_fields)     ! number of thresholds for each field
        integer istart_a(n_fields)   ! start time for each field              
        integer ipersist_a(n_fields) ! persistence flag for each field              
        logical l_persist
        character*2 c2_region
        character*10 c_thr

!       Specify what is being verified
        data var_a      
     1     /'SWI','TSF','DSF','USF','VSF','SSF','WSF'
     1     ,'T3' ,'W3' ,'TPW','R01','RTO','S8A'/     
        data nthr_a     /1,1,1,1,1,1,1,1,1,1,1,1,1/        
        data istart_a   /0,0,0,0,0,0,0,0,0,0,1,1,1/        
        data ipersist_a /0,0,0,0,0,0,0,0,0,1,0,0,1/        

        integer contable(0:1,0:1)

        integer maxthr
        parameter (maxthr=5)

      ! Never being used! Yuanfu covers this for now
        ! real cont_4d(NX_L,NY_L,NZ_L,maxthr)

        real obs_mean(maxbgmodels,0:max_fcst_times,max_regions)
        real fcst_mean(maxbgmodels,0:max_fcst_times,max_regions)
        real rms(maxbgmodels,0:max_fcst_times,max_regions)

        real obs_sum(maxbgmodels,0:max_fcst_times,max_regions)
        real fcst_sum(maxbgmodels,0:max_fcst_times,max_regions)
        real sumsq(maxbgmodels,0:max_fcst_times,max_regions)

        real obs_mean_comp(maxbgmodels,0:max_fcst_times,max_regions)
        real fcst_mean_comp(maxbgmodels,0:max_fcst_times,max_regions)
        real rms_comp(maxbgmodels,0:max_fcst_times,max_regions)

        real 
     1  frac_coverage(maxbgmodels,0:max_fcst_times,max_regions)
        real 
     1  frac_cvr_comp(maxbgmodels,0:max_fcst_times,max_regions)
        integer n(maxbgmodels,0:max_fcst_times,max_regions) 
        integer n_sum(maxbgmodels,0:max_fcst_times,max_regions) 

        integer nsuccess_m(maxbgmodels)
        integer nincomplete_m(maxbgmodels)
        integer incomplete_run_m(maxbgmodels)
        integer n_plot_times_m(maxbgmodels,0:max_fcst_times,n_fields)       
        integer i_plot_times_m(maxbgmodels,0:max_fcst_times,n_fields)       

       write(6,*)
       write(6,*)' Start subroutine verif_fcst_pt_composite...'
 
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

        rmiss = -99.9
        imiss = -999

        thresh_var = 20. ! lowest threshold for this variable

        write(6,*)' Time is ',i4time_sys,a9time
        write(6,*)' n_plot_times ',n_plot_times

        n_plot_times_m(:,:,:) = 1
        i_plot_times_m(:,:,:) = 0

        i4_initial = i4time_sys

!       Determine verification timing
        model_verif_intvl = max(laps_cycle_time,model_fcst_intvl)
        n_fcst_times = (model_fcst_len*60) / model_verif_intvl

        write(6,*)' model_verif_intvl = ',model_verif_intvl
        write(6,*)' n_fcst_times = ',n_fcst_times

        lun_in = 21

        lun_mem = 25

        lun_stats_in = 31
        lun_stats_out = 41
        lun_summary_out = 43

        do ifield = 1,n_fields

         if(ipersist_a(ifield) .eq. 1)then
          l_persist = .true.
         else
          l_persist = .false.
         endif

!        Get fdda_model_source and 'n_fdda_models' from static file
         call get_fdda_model_source(c_fdda_mdl_src,n_fdda_models
     1                                            ,istatus)

         if(l_persist .eqv. .true.)then
            n_fdda_models = n_fdda_models + 1
            c_fdda_mdl_src(n_fdda_models) = 'persistence'
            write(6,*)' Adding persistence to fdda_models'
         endif

         write(6,*)' n_fdda_models = ',n_fdda_models
         write(6,*)' c_fdda_mdl_src = '
     1             ,(c_fdda_mdl_src(m),m=1,n_fdda_models)

!        Update array of n_plot_times_m for available model exceptions
         do imodel = 2,n_fdda_models
             if(var_a(ifield) .eq. 'RTO')then
                 n_plot_times_m(imodel,0:max_fcst_times,ifield) = 0       

                 itime_exception  = 3600 / model_verif_intvl
                 if(itime_exception .le. max_fcst_times)then
                     n_plot_times_m(imodel,itime_exception,ifield) = 1
                 endif

                 itime_exception = 10800 / model_verif_intvl
                 if(itime_exception .le. max_fcst_times)then
                     n_plot_times_m(imodel,itime_exception,ifield) = 1
                 endif

                 itime_exception = 43200 / model_verif_intvl
                 if(itime_exception .le. max_fcst_times)then
                     n_plot_times_m(imodel,itime_exception,ifield) = 1
                 endif

                 itime_exception = 86400 / model_verif_intvl
                 if(itime_exception .le. max_fcst_times)then
                     n_plot_times_m(imodel,itime_exception,ifield) = 1
                 endif
             endif ! RTO field    

!            Accept just 3-hourly forecasts from the NAM model
             if(trim(c_fdda_mdl_src(imodel)) .eq. 'nam')then 
                 nam_fcst_intvl = 10800                           
                 do itime_fcst = 0,n_fcst_times
                     i4_fcst = itime_fcst*model_verif_intvl      
                     if(i4_fcst .ne. 
     1                 (i4_fcst/nam_fcst_intvl)*nam_fcst_intvl)then
                         n_plot_times_m(imodel,itime_fcst,ifield) = 0
                     endif
                 enddo
             endif ! NAM model

!            Accept just 6-hourly forecasts from the NAM-NH model
             if(trim(c_fdda_mdl_src(imodel)) .eq. 'nam-nh')then 
                 nam_fcst_intvl = 21600                           
                 do itime_fcst = 0,n_fcst_times
                     i4_fcst = itime_fcst*model_verif_intvl      
                     if(i4_fcst .ne. 
     1                 (i4_fcst/nam_fcst_intvl)*nam_fcst_intvl)then
                         n_plot_times_m(imodel,itime_fcst,ifield) = 0
                     endif
                 enddo
             endif ! NAM model

!            Accept just 1-hourly 3D/FUA fcsts from the HRRR model
             if(trim(c_fdda_mdl_src(imodel)) .eq. 'wrf-hrrr' .AND.
     1          (var_a(ifield) .eq. 'TPW' .or. 
     1           var_a(ifield) .eq. 'T3'  .or. 
     1           var_a(ifield) .eq. 'W3'       ) 
     1                                                             )then
                 mdl_fcst_intvl = 3600                            
                 do itime_fcst = 0,n_fcst_times
                     i4_fcst = itime_fcst*model_verif_intvl      
                     if(i4_fcst .ne. 
     1                 (i4_fcst/mdl_fcst_intvl)*mdl_fcst_intvl)then
                         n_plot_times_m(imodel,itime_fcst,ifield) = 0
                     endif
                 enddo
             endif ! NAM model

         enddo ! imodel

         call get_directory('verif',verif_dir,len_verif)

!        Initialize arrays
         obs_mean_comp = rmiss
         fcst_mean_comp = rmiss
         rms_comp = rmiss
         n_sum = 0
         obs_sum = 0.
         fcst_sum = 0.
         sumsq = 0.

         var_2d = var_a(ifield)
         call s_len(var_2d,lenvar)

         write(6,*)                                            
         write(6,*)' Processing field for composite stats ',var_2d               

         if(var_2d(1:lenvar) .eq. 'SWI')then
          model_cycle_time = 86400
         else
          model_cycle_time = model_cycle_time_in
         endif

         n_init_times = ((ndays * 86400) / model_cycle_time) - 1

         write(6,*)' ndays / n_init_times = ',ndays,n_init_times

         iregion = 1

         nmissing = 0
         nsuccess = 0
         nincomplete = 0

         nsuccess_m = 0
         nincomplete_m = 0

         frac_thr = 0.20
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

          plot_dir = verif_dir(1:len_verif)//var_2d(1:lenvar)
     1                                     //'/pt'
          len_plot = len_verif + 3 + lenvar ! + len_model

!         Read individual bias files
          write(6,*)' nthr before assigning = ',nthr
          nthr = nthr_a(ifield) ! nthr may be unset or have earlier been stepped on
          do idbz = 1,nthr

           rdbz = float(idbz*10) + 10
           write(c_thr,901)nint(rdbz)
 901       format(i2)

!          write GNUplot file for the time series of this model (region 1)
           stats_file_in = plot_dir(1:len_plot)//'/'
     1                                     //a9time_initial     
     1                                     //'.stats'     

           write(6,*)'stats_file_in = ',stats_file_in

           members_file  = verif_dir(1:len_verif)//'/'
     1                                     //'members.txt'      

           write(6,*)'members_file = ',members_file

           inquire(file=stats_file_in,exist=l_exist)
           if(.not. l_exist)then
               nmissing = nmissing + 1
               if(nmissing .le. nmissing_thr)then
                   write(6,*)' WARNING: file does not exist:'
     1                                                    ,stats_file_in
                   goto960
               else
                   write(6,*)' WARNING: file does not exist:'
     1                                                    ,stats_file_in
                   write(6,*)
     1  ' Skipping this field, too many missing initialization times...'       
     1                       ,nmissing,nmissing_thr                  
                   goto980
               endif
           endif ! l_exist

!          inquire(file=stats_file_in,exist=l_exist)
!          if(.not. l_exist)then
!              write(6,*)' ERROR: file does not exist:',stats_file_in       
!              goto980
!          endif ! l_exist

           open(lun_stats_in,file=stats_file_in,status='old')

!          Read comment with model member names
           read(lun_stats_in,51,err=958,end=958) cline
51         format(a)
           write(6,*)'cline = ',cline

           char = ' '
           call csplit(cline,c_fdda_mdl_hdr,nelems,n_fdda_models,char
     1                ,istatus)
           if(istatus .ne. 1)then
               write(6,*)' Read error in stats file header line...'
               goto958
           endif
           write(6,*)'c_fdda_mdl_hdr nelems/n_fdda_models = '
     1                              ,nelems,n_fdda_models

902        format('# ',30(1x,a))  

           if(l_col)then
!              Read mean and rms values
               do imodel = 2,n_fdda_models
                 if(imodel .gt. 2)then
                     read(lun_stats_in,*)                       
                     read(lun_stats_in,*)                       
                 endif

                 do itime_fcst = 0,n_fcst_times
                   i4_valid = i4_initial + itime_fcst*model_verif_intvl      
                   call cv_i4tim_asc_lp(i4_valid,a24time_valid_expected
     1                                 ,istatus)

                   read(lun_stats_in,911,err=912,end=912)
     1                  a24time_valid_file,    
     1                  obs_mean(imodel,itime_fcst,iregion),
     1                  fcst_mean(imodel,itime_fcst,iregion),
     1                  rms(imodel,itime_fcst,iregion)
911                format(a24,1x,3f10.3)
                   goto 913
912                write(6,*)' Read error in stats file...'
                   goto 958
913                continue

                   call left_justify(a24time_valid_file)
                   call left_justify(a24time_valid_expected)        

                   if(a24time_valid_file .ne. 
     1                a24time_valid_expected)then
                       write(6,*)
     1                 ' WARNING: imodel / a24time (expected/file)'
     1                 ,imodel,itime,a24time_valid_expected
     1                              ,a24time_valid_file     
                       goto958
                   endif

                 enddo ! itime_fcst

               enddo ! imodel

!              Test for missing data in all times/models                    
               do itime_fcst = istart_a(ifield),n_fcst_times

                 do imodel = 2,n_fdda_models

                   if(obs_mean(imodel,itime_fcst,iregion).eq.rmiss)then
                     n(imodel,itime_fcst,iregion) = 0     
                     i_good_timestep_model = 0
                   elseif(obs_mean(imodel,itime_fcst,iregion).eq.0.
     1                               .AND. 
     1                    var_2d .ne. 'R01' .and. var_2d .eq. 'RTO')then
                     write(6,914)imodel,itime_fcst,n_fcst_times       
914                  format(' ERROR: zero obs_mean'
     1                     ,' imodel/itime = ',2i4
     1                     ,' check # fcst times in file is ',i4)
                     write(6,*)' Setting n array to zero'
                     n(:,:,:) = 0     
                     i_good_timestep_model = 0
                   else ! valid data
                     n(imodel,itime_fcst,iregion) = 1 ! nobs
                     i_good_timestep_model = 1
                     i_plot_times_m(imodel,itime_fcst,ifield) = 1       
                   endif

!                  Flag as incomplete if missing during expected time     
                   if(i_good_timestep_model .eq. 0 .AND.
     1                n_plot_times_m(imodel,itime_fcst,ifield) .eq. 1 
     1                                                   )then
                     if(incomplete_run .eq. 0     .AND. 
     1                  itime_fcst .le. n_plot_times)then
                       write(6,916)init,a9time_initial,itime_fcst,imodel
916                    format(
     1                 ' WARNING: overall missing for init/time/model  '
     1                        ,i6,1x,a9,1x,i6,i6)                      
                       incomplete_run = 1 
                     endif

                     if(incomplete_run_m(imodel) .eq. 0     .AND. 
     1                  itime_fcst .le. n_plot_times)then
                       write(6,917)init,a9time_initial,itime_fcst,imodel
917                    format(
     1                 ' WARNING: missing N values for init/time/model '
     1                        ,i6,1x,a9,1x,i6,i6)                      
                       incomplete_run_m(imodel) = 1
                     endif

                   endif ! i_good_timestep_model = 0 .AND.
                         ! expect output from model at this time

                 enddo ! imodel
               enddo ! itime_fcst

!              Read members.txt file
               open(lun_mem,file=members_file,status='old')
               do imodel = 2,n_fdda_models
                 if(trim(c_model) .ne. 'persistence')then
                   if(c_model(1:3) .ne. 'lga')then
                       read(lun_mem,*)c_model               
                   endif
                   write(6,*)' c_model check (members/namelist): '
     1                                         ,imodel,c_model
     1                                         ,c_fdda_mdl_src(imodel)
!    1                                         ,c_fdda_mdl_hdr(imodel)

!                  Compare members.txt file and namelist fdda parms
                   if(trim(c_model) .ne.
     1                              trim(c_fdda_mdl_src(imodel)))then
                       write(6,*)
     1             ' ERROR: Models did not match (members.txt/namelist)'
                       close(lun_mem)
                       return
                   endif

!                  Compare rms file header and namelist fdda parms
                   if(trim(c_model) .ne.
     1                              trim(c_fdda_mdl_hdr(imodel)))then
                       write(6,*)
     1       ' WARNING: Models did not match (rms file header/namelist)'
                       close(lun_mem)
                       goto 958
                   endif

                 endif ! other than persistence from members.txt

               enddo ! imodel
               close(lun_mem)

           endif

940        close(lun_stats_in)

          enddo ! idbz

!         Write to log for informational purposes
          do imodel=2,n_fdda_models
             do itime_fcst = 0,n_fcst_times
                 write(6,950)c_fdda_mdl_src(imodel)
     1                    ,itime_fcst
     1                    ,n(imodel,itime_fcst,1)
     1                    ,n_sum(imodel,itime_fcst,1)                
     1                    ,obs_mean(imodel,itime_fcst,1)                
     1                    ,obs_sum(imodel,itime_fcst,1)                
     1                    ,fcst_mean(imodel,itime_fcst,1)                
     1                    ,fcst_sum(imodel,itime_fcst,1)                
     1                    ,rms(imodel,itime_fcst,1)                
     1                    ,sumsq(imodel,itime_fcst,1)                
950              format(' Obs counts/sums: ',a10,' itime_fcst' 
     1                                      ,i3,3x,2i6,3(f9.3,f9.1))
                 if(obs_mean(imodel,itime_fcst,1) .eq. 0.    
     1                               .AND. 
     1              var_2d .ne. 'R01' .and. var_2d .eq. 'RTO')then
                     write(6,*)' WARNING: obs_mean = 0. itime = '
     1                        ,itime_fcst
                 endif
             enddo ! itime_fcst
          enddo ! imodel

          nsuccess_m(:) = nsuccess_m(:) + (1 - incomplete_run_m(:))

          nincomplete = nincomplete + incomplete_run
          nincomplete_m(:) = nincomplete_m(:) + incomplete_run_m(:)

          nincomplete_t = 0
          do imodel=2,n_fdda_models
            if(trim(c_fdda_mdl_src(imodel)) .ne. 'advection')then
              nincomplete_t = nincomplete_t + incomplete_run_m(imodel)
            endif
          enddo 

          if(nincomplete_t .eq. 0 .OR. 
     1                             (l_req_all_mdls .eqv. .false.) )then
              write(6,*)' Accumulating sums for this run '
              nsuccess = nsuccess + 1

              where (n(:,:,:) .ne. imiss .AND. n(:,:,:) .gt. 0
     1                       .AND. obs_mean(:,:,:) .ne. rmiss
     1                       .AND. fcst_mean(:,:,:) .ne. rmiss
     1                                                           )                
               n_sum(:,:,:) = n_sum(:,:,:) + n(:,:,:)
               obs_sum(:,:,:) = 
     1         obs_sum(:,:,:)   + n(:,:,:)*obs_mean(:,:,:)
               fcst_sum(:,:,:) = 
     1         fcst_sum(:,:,:) + n(:,:,:)*fcst_mean(:,:,:)
               sumsq(:,:,:) = sumsq(:,:,:) + n(:,:,:)*(rms(:,:,:)**2)
              end where

          else
              write(6,*)' Not accumulating sums for this run'

          endif ! accumulate sums for this run

955       close(lun_stats_in)                                      

          write(6,956)a9time_initial,nincomplete_t,
     1                (incomplete_run_m(imodel),imodel=2,n_fdda_models)
956       format(' incomplete_run_m at ',a9,' is ',i3,4x,20i3)   

          goto 960 ! success for this time

!         Error Condition for this time
958       nmissing = nmissing + 1

960      enddo                    ! init (initialization time)

         write(6,*)'init neg ',n_sum(2,0,1)

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

!        Define and write to summary*.txt file
         summary_file_out = verif_dir(1:len_verif)//var_2d(1:lenvar)
     1                              //'/pt'
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
         close(lun_summary_out)

         if(nsuccess .lt. nsuccess_thr)then
             write(6,*)' Insufficient successful times to plot'
             goto 980
         endif

!        Calculate composite mean rms
         if(.true.)then
          do idbz = 1,nthr
           do imodel=2,n_fdda_models

             write(6,971)imodel,c_fdda_mdl_src(imodel)         
971          format(/' Writing sums for model ',i3,1x,a10)       

             do itime_fcst = 0,n_fcst_times

                 write(6,972)
     1                 n_sum(imodel,itime_fcst,iregion)
     1                ,obs_sum(imodel,itime_fcst,iregion)
     1                ,fcst_sum(imodel,itime_fcst,iregion)
     1                ,sumsq(imodel,itime_fcst,iregion)
     1                ,imodel,itime_fcst,iregion           
972              format('sums ',i8,3f9.1,3i4)

             enddo ! itime_fcst
           enddo ! imodel
          enddo ! idbz
         endif ! .true. / .false.

         where (n_sum(:,:,:) .ne. imiss .AND. n_sum(:,:,:) .gt. 0)
           obs_mean_comp(:,:,:) = obs_sum(:,:,:) / n_sum(:,:,:)
           fcst_mean_comp(:,:,:) = fcst_sum(:,:,:) / n_sum(:,:,:)
           rms_comp(:,:,:) = sqrt(sumsq(:,:,:) / n_sum(:,:,:))
         end where

         write(6,*)
         write(6,*)
     1    ' Write output composite obs/fcst/rms files ############# '                                           
         write(6,*)' nthr before assigning = ',nthr
         nthr = nthr_a(ifield) ! nthr may be unset or have earlier been stepped on
         do idbz = 1,nthr

           rdbz = float(idbz*10) + 10
           write(c_thr,901)nint(rdbz)

           i4_initial = i4time_sys                                   
           call make_fnam_lp(i4_initial,a9time_initial,istatus)

!          write GNUplot file for the time series of this model (region 1)
           stats_file_out = plot_dir(1:len_plot)//'/'
     1                                     //trim(compdir)//'/'
     1                                     //a9time_initial     
     1                                     //'.stats'     

           write(6,*)'stats_file_out = ',stats_file_out

           open(lun_stats_out,file=stats_file_out,status='unknown')

!          Write comment with model member names
           write(lun_stats_out,902)(trim(c_fdda_mdl_src(imodel))
     1                              ,imodel=2,n_fdda_models)  
           if(l_col)then

             do imodel=2,n_fdda_models

               if(imodel .gt. 2)then
                   write(6,*)                       
                   write(6,*)                       

                   write(lun_stats_out,*)                       
                   write(lun_stats_out,*)                       
               endif

!              Write mean and rms values
               do itime_fcst = 0,n_fcst_times
                   i4_valid = i4_initial + itime_fcst*model_verif_intvl      
                   call cv_i4tim_asc_lp(i4_valid,a24time_valid
     1                                 ,istatus)

                   if(l_time_outcoord_hhmm .eqv. .true.)then
                     i4_fcst = itime_fcst*model_verif_intvl      
                     i4_fcst_hh = i4_fcst / 3600
                     i4_fcst_mm = i4_fcst/60 - (i4_fcst_hh*60)

                     write(6,973)i4_fcst_hh,i4_fcst_mm,    
     1                 obs_mean_comp(imodel,itime_fcst,iregion),
     1                 fcst_mean_comp(imodel,itime_fcst,iregion), 
     1                 rms_comp(imodel,itime_fcst,iregion) 
973                  format(1x,i2.2,1x,i2.2,19x,3f10.3)

                     write(lun_stats_out,973)i4_fcst_hh,i4_fcst_mm,    
     1                 obs_mean_comp(imodel,itime_fcst,iregion),
     1                 fcst_mean_comp(imodel,itime_fcst,iregion), 
     1                 rms_comp(imodel,itime_fcst,iregion) 

                   else ! use ascii time format

                     write(6,974)a24time_valid,    
     1                 obs_mean_comp(imodel,itime_fcst,iregion),
     1                 fcst_mean_comp(imodel,itime_fcst,iregion), 
     1                 rms_comp(imodel,itime_fcst,iregion) 
974                  format(1x,a24,3f10.3)

                     write(lun_stats_out,974)a24time_valid,    
     1                 obs_mean_comp(imodel,itime_fcst,iregion),
     1                 fcst_mean_comp(imodel,itime_fcst,iregion), 
     1                 rms_comp(imodel,itime_fcst,iregion) 

                   endif ! l_time_outcoord_hhmm

               enddo ! itime_fcst

             enddo ! imodel

           endif

           close(lun_stats_out)

         enddo ! idbz

 980    enddo                     ! fields

       enddo ! i_period

 999   write(6,*)' End of subroutine verif_fcst_composite'

       return

       end

