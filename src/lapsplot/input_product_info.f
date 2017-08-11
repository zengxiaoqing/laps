
        subroutine input_product_info(
     1                              i4time_ref              ! I
     1                             ,laps_cycle_time         ! I
     1                             ,ndim                    ! I
     1                             ,c_prodtype              ! O
     1                             ,ext                     ! O
     1                             ,directory               ! O
     1                             ,a9time                  ! O
     1                             ,fcst_hhmm               ! O
     1                             ,i4_initial              ! O
     1                             ,i4_valid                ! O
     1                             ,istatus)                ! O

        include 'lapsparms.for' ! maxbgmodels

!       integer       maxbgmodels
!       parameter     (maxbgmodels=10)

        integer       n_fdda_models
        integer       l,len_dir,lfdda
        integer       istatus
        character*9   c_fdda_mdl_src(maxbgmodels)
        character*(*) directory
        character*(*) ext
        character*40  c_model
        character*10  cmds
        character*1   cansw
        character*150 c_filenames(1000)

        character*1 c_prodtype
        character*5 fcst_hhmm
        character*9 a9time

        logical l_parse

        write(6,*)' Subroutine input_product_info...'

        write(6,1)
 1      format('  Product type: analysis [a], bkgnd [b]'
     1        ,', balance [n], fcst [f], quit [q] ? ',$)        

        read(5,2)c_prodtype
 2      format(a)

        call upcase(c_prodtype,c_prodtype)

        if(c_prodtype .eq. 'A' .or. c_prodtype .eq. 'N')then
            i4_valid = i4time_ref
            call make_fnam_lp(i4_valid,a9time,istatus)
            write(6,*)' return from input_product_info ',a9time
            return

        elseif(c_prodtype .eq. 'Q')then
            write(6,*)' Quitting'
            return

        else
            if(c_prodtype .eq. 'B')then
                if(ndim .eq. 2)then
                    ext = 'lgb'
                else
                    ext = 'lga'
                endif

            elseif(c_prodtype .eq. 'F')then
                if(ndim .eq. 2)then
                    ext = 'fsf'
                else
                    ext = 'fua'
                endif

            endif

            call input_background_info(
     1                              ext                     ! I
     1                             ,directory,c_model       ! O
     1                             ,i4time_ref              ! I
     1                             ,laps_cycle_time         ! I
     1                             ,a9time                  ! O
     1                             ,fcst_hhmm               ! O
     1                             ,i4_initial              ! O
     1                             ,i4_valid                ! O
     1                             ,istatus)                ! O

        endif

        return
        end


        subroutine input_background_info(
     1                              ext                     ! I
     1                             ,directory,c_model       ! O
     1                             ,i4time_ref              ! I
     1                             ,laps_cycle_time         ! I
     1                             ,asc9_tim              ! O
     1                             ,fcst_hhmm               ! O
     1                             ,i4_initial              ! O
     1                             ,i4_valid                ! O
     1                             ,istatus)                ! O

        include 'lapsparms.for' ! maxbgmodels, max_background_files

!       integer       maxbgmodels
!       parameter     (maxbgmodels=10)

        integer       maxfiles
        parameter     (maxfiles=max_background_files)

        integer       n_fdda_models
        integer       l,len_dir,lfdda
        integer       istatus
        character*9   c_fdda_mdl_src(maxbgmodels)
        character*(*) directory
        character*(*) ext
        character*40  c_model
        character*10  cmds
        character*1   cansw
        character*150 c_filenames(maxfiles)

        character*5 fcst_hhmm
        character*9 asc9_tim, a9time
        character*14 a14_time

        logical l_parse

        write(6,*)' Subroutine input_background_info...'

        write(6,*)' Using ',ext(1:3),' file'

        istatus = 0

        call get_directory(ext,directory,len_dir)

        if(l_parse(ext,'lga') .or. l_parse(ext,'lgb'))then ! use LGA/LGB
            c_model = ' '
            write(6,*)' lga/lgb ext - setting c_model to blank'
            go to 900 
        endif

!       Get fdda_model_source from static file
        call get_fdda_model_source(c_fdda_mdl_src,n_fdda_models,istatus)

!       call get_file_names(directory,n_fdda_models,c_fdda_mdl_src
!    1                                             ,maxbgmodels,istatus)       

        n_fdda_models = n_fdda_models + 1
        c_fdda_mdl_src(n_fdda_models) = 'mm5'
        n_fdda_models = n_fdda_models + 1
        c_fdda_mdl_src(n_fdda_models) = 'mm5hot'

        call s_len(directory,len_dir)
        cansw='n'
        l=1

        if(n_fdda_models.eq.0)then
!          print*,'fdda is not turned on in static file'
!          return 
           write(6,*)' Assuming LGA only since n_fdda_models was zero'       
           n_fdda_models = 1
           c_fdda_mdl_src(1) = 'lga'
        endif

        write(6,*)' Available models are...'

        do l = 1,n_fdda_models
            call s_len(c_fdda_mdl_src(l),lfdda)
            if(c_fdda_mdl_src(l)(1:lfdda) .ne. 'lga')then
                write(6,*)' ',c_fdda_mdl_src(l)(1:lfdda)
            endif
        enddo ! l

        call s_len(c_fdda_mdl_src(1),lfdda)
        write(6,205)c_fdda_mdl_src(1)(1:lfdda),ext(1:3)
 205    format(/'  Enter model [e.g. ',a,'] for ',a3,' file: ',$)

        read(5,206)c_model
 206    format(a)

        write(6,*)' c_model = ',c_model

        call s_len(c_model,len_model)

        DIRECTORY=directory(1:len_dir)//c_model(1:len_model)//'/'

 900    continue

        call get_file_names(directory,nfiles,c_filenames
     1                     ,maxfiles,istatus)

        write(6,*)' Available files in ',trim(directory)
        if(nfiles .ge. 1)then
            do i = 1,nfiles
                call s_len(c_filenames(i),len_fname)
                write(6,*)c_filenames(i)(1:len_fname)
            enddo
        endif

        call       input_model_time(i4time_ref              ! I
     1                             ,laps_cycle_time         ! I
     1                             ,asc9_tim              ! O
     1                             ,fcst_hhmm               ! O
     1                             ,i4_initial              ! O
     1                             ,i4_valid                ! O
     1                                                            )

        istatus = 1

        if(i4_initial .eq. 0 .or. i4_valid .eq. 0)then ! find best fcst

            write(6,*)' Looking for best file'

            call get_best_fcst(maxfiles,i4time_ref,nfiles
     1                        ,c_filenames,i_best_file)
            if(i_best_file .gt. 0)then ! File for this ext exists with proper
               i = i_best_file
               call get_directory_length(c_filenames(i),lend)
               call get_time_length(c_filenames(i),lenf)
               a14_time = c_filenames(i)(lend+1:lenf)

               write(6,*)' Found file for: ',c_filenames(i)(lend+1:lenf)       
     1                                     ,' ',ext(1:6)

               call get_fcst_times(a14_time,i4_initial,i4_valid,i4_fn)
               write(6,*)' a14_time = ',a14_time
               call s_len(a14_time,length_fcst)
               fcst_hhmm = a14_time(10:length_fcst)
               call make_fnam_lp(i4_valid,asc9_tim,istatus)
               write(6,*)' Valid time = ',asc9_tim


               write(6,*)' Found file for: ',c_filenames(i)(lend+1:lenf)
     1                                           ,' ',ext(1:6)

            else
               call make_fnam_lp(i4time_ref,a9time,istatus)
               write(6,*)' Could not find best file valid at: '
     1                   ,a9time
               fcst_hhmm = 'N/A '
               istatus = 0

            endif

        endif

        write(6,*)' Exit input_background_info, ASCII valid time = ',
     1            asc9_tim

        return
        end

        subroutine input_model_time(i4time_ref              ! I
     1                             ,laps_cycle_time         ! I
     1                             ,asc9_tim              ! O
     1                             ,fcst_hhmm               ! O
     1                             ,i4_initial              ! O
     1                             ,i4_valid                ! O
     1                                                            )

        character*5 fcst_hhmm
        character*9 asc9_tim, a9time
        character*14 a14_time

 1200   write(6,211)
 211    format(/'  Enter yydddhhmm[H]HHMM or [H]HHMM for file, '
     1         ,'or blank for best fcst: ',$)

        read(5,221)a14_time
 221    format(a14)

        call s_len(a14_time,len_time)

        if(len_time .eq. 13 .or. len_time .eq. 14)then ! yydddhhmmHHMM/HHHMM
                write(6,*)' len_time = ',len_time
                call get_fcst_times(a14_time,i4_initial,i4_valid,i4_fn)
                write(6,*)' a14_time = ',a14_time
                if(len_time .eq. 13)then ! yydddhhmmHHMM
                    fcst_hhmm = a14_time(10:13)
                else                     ! yydddhhmmHHHMM
                    fcst_hhmm = a14_time(10:14)
                endif
                call make_fnam_lp(i4_valid,asc9_tim,istatus)
                write(6,*)' Valid time = ',asc9_tim

        elseif(len_time .eq. 4 .or. len_time .eq. 5)then
                write(6,*)' len_time = ',len_time

                i4time_plot = i4time_ref ! / laps_cycle_time 
!    1                                     * laps_cycle_time       
                call make_fnam_lp(i4time_plot,asc9_tim,istatus)
                write(6,*)' Valid time = ',asc9_tim

                fcst_hhmm = a14_time(1:len_time)

              ! Get fcst interval
                a14_time = asc9_tim//fcst_hhmm
                call get_fcst_times(a14_time,I4TIME,i4_valid,i4_fn) 
                i4_interval = i4_valid - I4TIME
                i4_initial = I4TIME - i4_interval ! Reset initial time
                i4_valid = i4_valid - i4_interval
                call make_fnam_lp(i4_initial,a9time,istatus)

                a14_time = a9time//fcst_hhmm
                write(6,*)' Modified a14_time = ',a14_time

        elseif(len_time .eq. 0)then
                write(6,*)' Input fcst time was blank'
                i4_initial = 0
                i4_valid = 0

        else
                write(6,*)' Try again, len_time = ',len_time
                goto1200

        endif

        return
        end

        subroutine input_level(lun,k_level,k_mb,pres_3d,NX_L,NY_L,NZ_L)       

!       Round to the nearest pressure level

        real pres_3d(NX_L,NY_L,NZ_L)

        character*40 vert_grid

        call get_vertical_grid(vert_grid,istatus)

        icen = NX_L/2
        jcen = NY_L/2

        read(lun,*)k_level
        k_mb = k_level
        if(k_level .gt. 0 .and. vert_grid .eq. 'PRESSURE')then
            pressure = float(k_level*100)
            k_level = nint(rlevel_of_field(pressure,pres_3d
     1                       ,NX_L,NY_L,NZ_L,icen,jcen,istatus))
            k_mb    = nint(pres_3d(icen,jcen,k_level) / 100.)
        endif

        return
        end
