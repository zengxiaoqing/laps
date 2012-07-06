          
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
        character*150 bias_file_in, ets_file_in
        character*150 bias_file_out, ets_file_out
        character*150 summary_file_out
        character*10 compdir

        integer n_fields
        parameter (n_fields=2)
        character*10 ext_anal_a(n_fields), ext_fcst_a(n_fields)
        character*10 var_a(n_fields)
        integer nthr_a(n_fields) ! number of thresholds for each field
        character*2 c2_region
        character*10 c_thr

!       Specify what is being verified
!       data ext_fcst_a /'fua'/ ! 3-D
!       data ext_anal_a /'lps'/ ! 3-D reflectivity
!       data var_a      /'REF'/ ! 3-D reflectivity       
!       data nthr_a     /5/        

        data ext_fcst_a /'fua','fsf'/ ! 3-D / composite ref
        data ext_anal_a /'lps','lmr'/ ! 3-D / composite ref
        data var_a      /'REF','LMR'/ ! 3-D / composite ref
        data nthr_a     /5,5/        

        integer contable(0:1,0:1)

        integer maxthr
        parameter (maxthr=5)

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
        integer 
     1  n_sum(maxbgmodels,0:max_fcst_times,max_regions,maxthr,0:1,0:1)

        integer nmissing_m(maxbgmodels)
        integer nsuccess_m(maxbgmodels)
        integer nincomplete_m(maxbgmodels)
        integer incomplete_run_m(maxbgmodels)

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

!       Initialize arrays
        bias = -999.
        ets = -999.
        frac_coverage = -999.
        n_sum = 0

        thresh_var = 20. ! lowest threshold for this variable

        write(6,*)' Time is ',i4time_sys,a9time
        write(6,*)' n_plot_times ',n_plot_times

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

        write(6,*)' ndays / n_init_times = ',ndays,n_init_times

        do ifield = 1,n_fields

         var_2d = var_a(ifield)
         call s_len(var_2d,lenvar)

         iregion = 1

         nmissing = 0
         nsuccess = 0
         nincomplete = 0

         nmissing_m = 0
         nsuccess_m = 0
         nincomplete_m = 0

         frac_thr = 0.15
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

           rdbz = float(idbz*10) + 10
           write(c_thr,901)nint(rdbz)
 901       format(i2)

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
               nmissing = nmissing + 1
               nmissing_m = nmissing_m + 1
!              if(nmissing .le. nmissing_thr)then
                   write(6,*)' WARNING: file does not exist:'
     1                                                     ,bias_file_in
                   goto960
!              else
!                  write(6,*)' ERROR: file does not exist:',bias_file_in       
!                  write(6,*)
!    1  ' Skipping this field, too many missing initialization times...'       
!    1                       ,nmissing_thr                  
!                  goto980
!              endif
           endif ! l_exist

           inquire(file=ets_file_in,exist=l_exist)
           if(.not. l_exist)then
               write(6,*)' ERROR: file does not exist:',ets_file_in       
               goto980
           endif ! l_exist

           open(lun_bias_in,file=bias_file_in,status='old')
           open(lun_ets_in,file=ets_file_in,status='old')

!          Read comment with model member names
           read(lun_bias_in,*)
           read(lun_ets_in,51) cline
 51        format(a)
           write(6,*)'cline = ',cline

           char = ' '
           call csplit(cline,c_fdda_mdl_hdr,nelems,n_fdda_models,char)
           write(6,*)'c_fdda_mdl_hdr nelems/n_fdda_models = '
     1                              ,nelems,n_fdda_models

902        format('# ',30(1x,a))  

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
     1                 ,' WARNING: imodel / a24time (expected/file-1)'
     1                 ,imodel,itime,a24time_valid_expected            
     1                              ,a24time_valid              
                       goto960
                   endif

                   read(lun_ets_in,911)a24time_valid,    
     1                 (ets(imodel,itime_fcst,iregion,idbz)
     1                              ,imodel=2,n_fdda_models)     
                   call left_justify(a24time_valid)
                   if(a24time_valid .ne.
     1                a24time_valid_expected)then
                       write(6,*)
     1                 ,' WARNING: imodel / a24time (expected/file-2)'
     1                 ,imodel,itime,a24time_valid_expected           
     1                              ,a24time_valid         
                       goto960
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
                   goto914
913                write(6,*)' read ERROR in bias file N vals '
                   write(6,*)' check if model config has changed'
                   write(6,*)in,jn,itime_fcst
                   write(6,*)cline
                   goto 955
914                write(6,915)in,jn,itime_fcst,
     1                 (n(imodel,itime_fcst,iregion,idbz,in,jn)
     1                              ,imodel=2,n_fdda_models)     
915                format('in,jn,itime,n',i2,i2,i4,20i10)

                   call left_justify(a24time_valid)
                   if(a24time_valid .ne.
     1                a24time_valid_expected)then
                       write(6,*)
     1                 ,' WARNING: imodel / a24time (expected/file-3)'
     1                 ,imodel,itime,a24time_valid_expected           
     1                              ,a24time_valid              
                       goto960
                   endif
                 enddo ! itime_fcst
               enddo ! jn
               enddo ! in

!              Test for missing data in all times/models for this dbz
               do itime_fcst = 0,n_fcst_times
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

!                  Flag as incomplete if missing during the first 12 hours
                   if(i_good_timestep_model .eq. 0)then
                     if(incomplete_run .eq. 0     .AND. 
     1                  itime_fcst .le. n_plot_times)then
                       write(6,916)init,itime_fcst,a9time_initial,imodel
916                    format(
     1                 ' WARNING: missing N values for init/time/model '
     1                        ,2i6,1x,a9,1x,i6)                      
                       incomplete_run = 1 
                     endif

                     if(incomplete_run_m(imodel) .eq. 0     .AND. 
     1                  itime_fcst .le. n_plot_times)then
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
     1                 ,' WARNING: imodel / a24time (expected/file-4)'
     1                 ,imodel,itime,a24time_valid_expected           
     1                              ,a24time_valid                           
                       goto960
                   endif
               enddo ! itime_fcst

!              Read members.txt file
               open(lun_mem,file=members_file,status='unknown')
               do imodel = 2,n_fdda_models
                   if(c_model(1:3) .ne. 'lga')then
                       read(lun_mem,*)c_model               
                   endif
                   write(6,*)' c_model check: ',c_model
     1                                         ,c_fdda_mdl_src(imodel)
     1                                         ,c_fdda_mdl_hdr(imodel)

!                  Compare members.txt file and namelist fdda parms
                   if(trim(c_model) .ne. 
     1                              trim(c_fdda_mdl_src(imodel)))then
                       write(6,*)' Models did not match (members.txt)'
                       return
                   endif

!                  Compare ets file header and namelist fdda parms
                   if(trim(c_model) .ne. 
     1                              trim(c_fdda_mdl_hdr(imodel)))then
                       write(6,*)' Models did not match (ets header)'
                       return
                   endif
               enddo ! imodel
               close(lun_mem)

           endif

940        close(lun_bias_in)
           close(lun_ets_in)

          enddo ! idbz 

          where (n(:,:,:,:,:,:) .ne. imiss)                
           n_sum(:,:,:,:,:,:) = n_sum(:,:,:,:,:,:) + n(:,:,:,:,:,:)
          end where

!         Write to log for informational purposes
          do imodel=2,n_fdda_models
             do itime_fcst = 0,n_fcst_times
                 write(6,950)c_fdda_mdl_src(imodel)
     1                    ,itime_fcst
     1                    ,n(imodel,itime_fcst,1,1,1,1)
     1                ,n_sum(imodel,itime_fcst,1,1,1,1)                
950              format(' Correct negs: ',a10,' itime_fcst' ,i3,3x,2i10)
             enddo ! itime_fcst
          enddo ! imodel

          nsuccess = nsuccess + 1
          nsuccess_m(:) = nsuccess_m(:) + (1 - incomplete_run_m(:))

          nincomplete = nincomplete + incomplete_run
          nincomplete_m(:) = nincomplete_m(:) + incomplete_run_m(:)

955       close(lun_bias_in)                                      
          close(lun_ets_in)                                         

          nincomplete_t = 0
          do imodel=2,n_fdda_models
              nincomplete_t = nincomplete_t + incomplete_run_m(imodel)
          enddo 

          write(6,956)a9time_initial,nincomplete_t,
     1                (incomplete_run_m(imodel),imodel=2,n_fdda_models)
956       format(' incomplete_run_m at ',a9,' is ',i3,4x,20i3)   

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

         if(nruns_plotted .ge. 2)then
             l_plot_criteria = .true.
         else
             l_plot_criteria = .false.
         endif

         write(6,*)'nruns_plotted / l_plot_criteria = ',nruns_plotted
     1                                                 ,l_plot_criteria

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
         close(lun_summary_out)

         if(nsuccess .lt. nsuccess_thr)then
             write(6,*)' Insufficient successful times to plot'
             goto 980
         endif

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
     1                  ,frac_obs                                      ! O
     1                  ,frac_fcst                                     ! O
     1                  ,frac_cvr_comp(imodel,itime_fcst,iregion,idbz) ! O
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

         write(6,*)
         write(6,*)
     1        ' Write output composite bias/ets files ############# '                                           
         write(6,*)' nthr before assigning = ',nthr
         nthr = nthr_a(ifield) ! nthr may be unset or have earlier been stepped on
         do idbz = 1,nthr

           rdbz = float(idbz*10) + 10
           write(c_thr,901)nint(rdbz)
 
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
           write(lun_bias_out,902)(trim(c_fdda_mdl_src(imodel))
     1                              ,imodel=2,n_fdda_models)  
           write(lun_ets_out,902)(trim(c_fdda_mdl_src(imodel))
     1                              ,imodel=2,n_fdda_models)  
           if(l_col)then
!              Write bias and ets values
               do itime_fcst = 0,n_fcst_times
                   i4_valid = i4_initial + itime_fcst*model_verif_intvl      
                   call cv_i4tim_asc_lp(i4_valid,a24time_valid
     1                                 ,istatus)

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
                   i4_valid = i4_initial + itime_fcst*model_verif_intvl      
                   call cv_i4tim_asc_lp(i4_valid,a24time_valid
     1                                 ,istatus)
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
                   i4_valid = i4_initial + itime_fcst*model_verif_intvl      
                   call cv_i4tim_asc_lp(i4_valid,a24time_valid
     1                                 ,istatus)
                   write(lun_bias_out,923)a24time_valid,    
     1                 (frac_cvr_comp                
     1                   (imodel,itime_fcst,iregion,idbz)
     1                              ,imodel=2,n_fdda_models)     
               enddo ! itime_fcst

           endif

           close(lun_bias_out)
           close(lun_ets_out)

         enddo ! idbz

 980    enddo                     ! fields

       enddo ! i_period

 999   write(6,*)' End of subroutine verif_radar_composite'

       return

       end


        subroutine csplit(line,carray,nelems,maxelems,char)

!       Routine to split a character array into segments

        character*(*) line
        character*1 char
        character*30 carray(maxelems)
        
        integer istart(maxelems)
        integer iend(maxelems)

        lenline = len(line)
        nelems = 0

        lenelem = len(carray(1))

!       Beginning of line might start an element
        if(line(1:1) .ne. char)then
            istart(1) = 1
        endif

        do i = 1,lenline-1

!           Check for start of string
            if(line(i:i) .eq. char .and. line(i+1:i+1) .ne. char)then
                if(nelems+1 .gt. maxelems)then
                    write(6,*)
     1              ' Error: nelems+1 > maxelems',nelems+1,maxelems,i
     1                     ,line(i:i+1)
                    write(6,*)line
                    stop
                endif
                istart(nelems+1)= i+1
            endif

!           Check for end of string
            if(line(i:i) .ne. char .and. line(i+1:i+1) .eq. char)then
                nelems = nelems + 1
                iend(nelems)= i+1
 
                if(iend(nelems) - istart(nelems) .ge. lenelem)then
                    write(6,*)' Error: element is too long',
     1                        istart(nelems),iend(nelems)
                    write(6,*)line
                    stop
                endif 

                carray(nelems) = line(istart(nelems):iend(nelems))
            endif

        enddo ! i

        write(6,*)' Elements in csplit are:'
        do i = 1,nelems
            write(6,*)i,' ',carray(i)
        enddo ! i

        return
        end
