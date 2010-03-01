
        program verif_radar_main

        character*9 a9time

        call get_systime(i4time,a9time,istatus)
        if(istatus .ne. 1)go to 999

        write(6,*)' systime = ',a9time

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
          
        call verif_radar(i4time,a9time,laps_cycle_time,
     1                  NX_L,NY_L,
     1                  NZ_L,
     1                  r_missing_data,
     1                  j_status)

999     continue

        end
          
        subroutine verif_radar(i4time_sys,a9time,laps_cycle_time,
     1                  NX_L,NY_L,
     1                  NZ_L,
     1                  r_missing_data,
     1                  j_status)

        real var_anal_3d(NX_L,NY_L,NZ_L)
        real var_fcst_3d(NX_L,NY_L,NZ_L)
        real rqc(NX_L,NY_L)

        logical lmask_and_3d(NX_L,NY_L,NZ_L)
        logical lmask_rqc_3d(NX_L,NY_L,NZ_L)

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
        character*9 a9time,a9time_valid
        character*150 hist_dir, cont_dir, verif_dir
        character*150 hist_file

        integer n_fields
        parameter (n_fields=1)
        character*10 ext_anal_a(n_fields), ext_fcst_a(n_fields)
        character*10 var_a(n_fields)
        integer nthr_a(n_fields) ! number of thresholds for each field
        character*2 c2_region

        data ext_fcst_a /'fua'/        
        data ext_anal_a /'lps'/        
        data var_a      /'REF'/        
        data nthr_a     /5/        

        integer contable(0:1,0:1)

        integer maxthr
        parameter (maxthr=5)

        real cont_4d(NX_L,NY_L,NZ_L,maxthr)

        thresh_var = 20. ! lowest threshold for this variable

        i4_initial = i4time_sys

        lun_in = 21

!       Get fdda_model_source from static file
        call get_fdda_model_source(c_fdda_mdl_src,n_fdda_models,istatus)

        write(6,*)' n_fdda_models = ',n_fdda_models
        write(6,*)' c_fdda_mdl_src = '
     1            ,(c_fdda_mdl_src(m),m=1,n_fdda_models)

!       Read in data file with region points
        call read_region_info(maxbgmodels,max_fcst_times,max_regions
     1                       ,n_models,n_fcst_times,n_regions
     1                       ,il,ih,jl,jh,lun_in)

        if(n_fdda_models .ne. n_models + 1)then
            write(6,*)' ERROR n_models differs from n_fdda_models '
     1                       ,n_models,n_fdda_models
            stop
        endif

        do ifield = 1,n_fields

         var_2d = var_a(ifield)
         call s_len(var_2d,lenvar)

         do imodel = 1,n_fdda_models

          c_model = c_fdda_mdl_src(imodel)

          if(c_model(1:3) .ne. 'lga')then

            write(6,*)' Processing model ',c_model

            call s_len(c_model,len_model)

            call get_directory('verif',verif_dir,len_verif)

            hist_dir = verif_dir(1:len_verif)//var_2d(1:lenvar)
     1                                       //'/hist/'
     1                                       //c_model(1:len_model)
            len_hist = len_verif + 6 + lenvar + len_model

            cont_dir = verif_dir(1:len_verif)//var_2d(1:lenvar)
     1                                       //'/cont/'
     1                                       //c_model(1:len_model)
     1                                       //'/'
!           len_cont = len_verif + 6 + lenvar + len_model

            do itime_fcst = 0,n_fcst_times

              itime = itime_fcst + 1

              do iregion = 1,n_regions ! 1 for testing

                ilow  = il(imodel,itime,iregion)
                ihigh = ih(imodel,itime,iregion)
                jlow  = jl(imodel,itime,iregion)
                jhigh = jh(imodel,itime,iregion)

                i4_valid = i4_initial + itime_fcst * laps_cycle_time

                call make_fnam_lp(i4_valid,a9time_valid,istatus)

                write(6,*)
                write(6,*)' Histograms for forecast time step '
     1                   ,itime_fcst,' Region = ',iregion

                lun_out = 11

                write(c2_region,1)iregion
 1              format(i2.2)

!               Add c_model to this?
                hist_file = hist_dir(1:len_hist)//'/'//a9time_valid       
     1                                          //'_'//c2_region     
     1                                          //'.hist'     

                write(6,*)'hist_file = ',hist_file

                open(lun_out,file=hist_file,status='unknown')

                if(iregion .eq. 1)then

!                 Read analyzed reflectivity
                  ext = ext_anal_a(ifield)
                  call get_laps_3d(i4_valid,NX_L,NY_L,NZ_L
     1             ,ext,var_2d,units_2d,comment_2d,var_anal_3d,istatus)
                  if(istatus .ne. 1)then
                        write(6,*)' Error reading 3D REF Analysis'
                        return
                  endif

                  if(var_2d .eq. 'REF')then ! also read radar quality
                      ext = 'lcv'
                      call get_laps_2d(i4_valid,ext,'RQC',units_2d
     1                                ,comment_2d,NX_L,NY_L,rqc,istatus)
                      if(istatus .ne. 1)then
                          write(6,*)' Error reading 2D RQC Analysis'
                          return
                      endif
                  endif

!                 Read forecast reflectivity
                  ext = ext_fcst_a(ifield)
                  call get_directory(ext,directory,len_dir)
                  DIRECTORY=directory(1:len_dir)//c_model(1:len_model)
     1                                          //'/'

                  call get_lapsdata_3d(i4_initial,i4_valid
     1                          ,NX_L,NY_L,NZ_L       
     1                          ,directory,var_2d
     1                          ,units_2d,comment_2d,var_fcst_3d
     1                          ,istatus)
                  if(istatus .ne. 1)then
                       write(6,*)' Error reading 3D REF Forecast'
                       return
                  endif

!                 Calculate "and" mask
                  n_rqc_and = 0
                  n_rqc_and_pot = 0
                  do k = 1,NZ_L
                  do i = 1,NX_L
                  do j = 1,NY_L
                     lmask_and_3d(i,j,k) = .false.
                     if(var_anal_3d(i,j,k) .ne. r_missing_data .and. 
     1                  var_anal_3d(i,j,k) .ge. thresh_var .and.
     1                  var_fcst_3d(i,j,k) .ne. r_missing_data .and.
     1                  var_fcst_3d(i,j,k) .ge. thresh_var        )then
                         lmask_and_3d(i,j,k) = .true.
                     endif
                     if(lmask_and_3d(i,j,k))then
                         n_rqc_and_pot = n_rqc_and_pot + 1
                         if(rqc(i,j) .ne. 3.0)then
                             n_rqc_and = n_rqc_and + 1
                             lmask_and_3d(i,j,k) = .false.
                         endif
                     endif
                  enddo ! j
                  enddo ! i
                  enddo ! k

!                 Calculate "rqc" mask
                  lmask_rqc_3d = .true.
                  n_rqc_all = 0
                  n_rqc_all_pot = 0
                  do k = 1,NZ_L
                  do i = 1,NX_L
                  do j = 1,NY_L
                    n_rqc_all_pot = n_rqc_all_pot + 1
                    if(rqc(i,j) .ne. 3.0)then
                        n_rqc_all = n_rqc_all + 1
                        lmask_rqc_3d(i,j,k) = .false.
                    endif
                  enddo ! j
                  enddo ! i
                  enddo ! k

                  write(6,*)' Time, # of QC (2D) points (and/all)'
     1                     ,itime_fcst,n_rqc_and,n_rqc_all
                  write(6,*)' Time, % of QC (2D) points (and/all)'
     1                     ,itime_fcst
     1                     ,float(n_rqc_and)/float(n_rqc_and_pot) * 100.
     1                     ,float(n_rqc_all)/float(n_rqc_all_pot) * 100.

                  if(n_rqc_all .eq. n_rqc_all_pot)then
                      write(6,*)' WARNING, 100% failed 3D radar QC test'
                      write(6,*)' Radar data is apparently all 2D'
                      write(6,*)' Histograms will be zeroed out'
                  endif

                endif ! iregion .eq. 1

                write(lun_out,*)
                write(lun_out,*)' REGION/IRANGE/JRANGE ',iregion
     1                         ,ilow,ihigh,jlow,jhigh

                write(lun_out,*)
                write(lun_out,*)' NO mask is in place'

                write(lun_out,*)
                write(lun_out,*)
     1                ' Calling radarhist for analysis at ',a9time_valid 
                call radarhist(NX_L,NY_L,NZ_L,var_anal_3d
     1                        ,ilow,ihigh,jlow,jhigh
     1                        ,lmask_rqc_3d,lun_out)       

                write(lun_out,*)
                write(lun_out,*)' Calling radarhist for',itime_fcst
     1                    ,' hr forecast valid at ',a9time_valid
                call radarhist(NX_L,NY_L,NZ_L,var_fcst_3d
     1                        ,ilow,ihigh,jlow,jhigh
     1                        ,lmask_rqc_3d,lun_out)

                write(lun_out,*)
                write(lun_out,*)
     1    ' 3-D AND mask is in place with dbz threshold of ',thresh_var       

                write(lun_out,*)
                write(lun_out,*)
     1                ' Calling radarhist for analysis at ',a9time_valid       
                call radarhist(NX_L,NY_L,NZ_L,var_anal_3d
     1                        ,ilow,ihigh,jlow,jhigh
     1                        ,lmask_and_3d,lun_out)

                write(lun_out,*)
                write(lun_out,*)' Calling radarhist for',itime_fcst
     1               ,' hr forecast valid at ',a9time_valid
                call radarhist(NX_L,NY_L,NZ_L,var_fcst_3d
     1                        ,ilow,ihigh,jlow,jhigh
     1                        ,lmask_and_3d,lun_out)

                nthr = nthr_a(ifield)

!               Calculate contingency tables
!               Radar QC (rqc) should be added perhaps via l_mask_all
!                              (or call l_mask_all l_mask_qconly)
                do idbz = 1,nthr
                  rdbz = float(idbz*10) + 10

                  write(lun_out,*)
                  write(lun_out,*)' Calculate contingency table for '
     1                           ,rdbz,' dbz'
                  write(lun_out,*)' region = ',iregion
     1                           ,ilow,ihigh,jlow,jhigh
                  call contingency_table(var_anal_3d,var_fcst_3d     ! I
     1                                  ,NX_L,NY_L,NZ_L,rdbz,lun_out ! I
     1                                  ,ilow,ihigh,jlow,jhigh       ! I
     1                                  ,lmask_rqc_3d                ! I
     1                                  ,contable)                   ! O


!                 Calculate/Write Skill Scores
                  call skill_scores(contable,lun_out)

                  if(iregion .eq. 1)then

!                     Calculate Contingency Table (3-D)
                      call calc_contable_3d(i4_initial,i4_valid
     1                       ,var_anal_3d,var_fcst_3d
     1                       ,rdbz,NX_L,NY_L,NZ_L                    ! I
     1                       ,lmask_rqc_3d,r_missing_data            ! I
     1                       ,cont_4d(1,1,1,idbz))                   ! O

                  endif ! iregion = 1

                enddo ! idbz

              enddo ! iregion

!             Write Contingency Tables (3-D)
              call put_contables(i4_initial,i4_valid,nthr
     1                        ,cont_4d,NX_L,NY_L,NZ_L,cont_dir)

              close (lun_out) 

            enddo ! itime_fcst

          endif ! c_model .ne. lga

         enddo ! model

        enddo ! fields

 999    write(6,*)' End of subroutine verif_radar'

        return

        end


        subroutine read_region_info(
     1                        maxbgmodels,max_fcst_times,max_regions
     1                       ,n_models,n_fcst_times,n_regions
     1                       ,il,ih,jl,jh,lun_in)

        integer il(maxbgmodels,max_fcst_times,max_regions)
        integer ih(maxbgmodels,max_fcst_times,max_regions)
        integer jl(maxbgmodels,max_fcst_times,max_regions)
        integer jh(maxbgmodels,max_fcst_times,max_regions)

        character*150 static_dir,static_file

!       Read in data file with region points
        call get_directory('static',static_dir,len_static)
        static_file = static_dir(1:len_static)//'/verif_regions.dat'
        open(lun_in,file=static_file,status='old')

        read(lun_in,*)n_regions
        read(lun_in,*)n_fcst_times
        read(lun_in,*)n_models

        do ir = 1,n_regions
            do if = 1,n_fcst_times
                read(lun_in,*)
                read(lun_in,*)i_fcst_time

                if(if .ne. (i_fcst_time+1))then
                    write(6,*)' ERROR in read_region_info '
                    write(6,*)' if differs from i_fcst_time+1: '
     1                        ,if,i_fcst_time,i_fcst_time+1
                    stop
                endif

                do im = 1,n_models
                    read(lun_in,*)il(im,if,ir),jl(im,if,ir)
                    read(lun_in,*)ih(im,if,ir),jh(im,if,ir)
                enddo ! im

            enddo ! if
        enddo ! ir
        close(lun_in)

        return
        end
