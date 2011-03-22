
        program ensemble_mean_main

! Program calculates ensmble mean as well as time mean
! Isidora Jankov, Linda Wharton and Steve Albers July 2010


        character*9 a9time
	character*300 laps_data_root,dir_t,filenamet
        character*250 stlaps
        integer      NX_L,NY_L,NZ_L,i4time,i,istatus 
        integer      model_cycle_time_sec, len_dir_t
        real         r_missing_data
        integer, parameter :: nprmax=150
	real, dimension(nprmax) :: pressures
        integer, parameter :: lun=120

	include 'lapsparms.for'

! Parameter for model cycle interval
        model_cycle_time_sec = 10800
!

	call getarg(1,laps_data_root)

!       call get_systime(i4time,a9time,istatus)
!       if(istatus .ne. 1)go to 999
!       write(6,*)' systime = ',a9time

        call get_directory('time',dir_t,istatus)
        call s_len(dir_t,len_dir_t)
        filenamet = dir_t(1:len_dir_t)//'/modeltime.dat'
        open(lun,file=filenamet,status='old')
        read(lun,*)a9time
        close(lun)
        call i4time_fname_lp(a9time,i4time,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'Error getting i4time ',a9time
           go to 999
        endif

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

	call get_pres_1d(i4time,NZ_L,pressures,istatus)
        if (istatus .ne. 1) then
           write (6,*) 'Error getting 1d pressures'
           go to 999
        endif

	

! Determine number of isobaric levels by checking pressure data.
!   (Assume there is at least one level, and that data is ordered correctly)

        do i=1,NZ_L
           if (pressures(i+1) > 200000.) exit
        enddo



        call ensemble_mean(i4time,model_cycle_time_sec,
     1                  NX_L,NY_L,
     1                  NZ_L,
     1                  pressures,	
     1                  r_missing_data,
     1                  istatus)

999     continue

        end

!******************************************************************
          
       subroutine ensemble_mean(i4time,model_cycle_time_sec,
     1                  NX_L,NY_L,
     1                  NZ_L,
     1                  pressures,	
     1                  r_missing_data,
     1                  istatus)


       include 'lapsparms.for'

       integer, allocatable, dimension(:) :: lvls3d
       integer, allocatable, dimension(:) :: lvls2d

       character(len=3),   allocatable, dimension(:) :: name3d
       character(len=10),  allocatable, dimension(:) :: units3d
       character(len=4),   allocatable, dimension(:) :: lvltype3d
       character(len=132), allocatable, dimension(:) :: com3d

       character(len=3),   allocatable, dimension(:) :: name2d
       character(len=10),  allocatable, dimension(:) :: units2d
       character(len=4),   allocatable, dimension(:) :: lvltype2d
       character(len=132), allocatable, dimension(:) :: com2d


       integer,intent(in) ::   NX_L,NY_L,NZ_L

        integer n_fields
        parameter (n_fields=4)

        real var_fcst_3d(NX_L,NY_L,NZ_L)
        real var_fcst_2d(NX_L,NY_L)
	real cdl_levels(NZ_L)
        real pressures(NZ_L) 

        real mean_fcst_3d(NX_L,NY_L,NZ_L)
        real mean_fcst_2d(NX_L,NY_L)

        real mean_fcst_3d_model(NX_L,NY_L,NZ_L*n_fields)
        real mean_fcst_2d_model(NX_L,NY_L,n_fields)

!        logical write_to_lapsdir

        integer       i4time, i4_initial, i4_valid
!        integer       maxbgmodels
!        parameter     (maxbgmodels=10)
        character*30  c_fdda_mdl_src(maxbgmodels)
        character*300 lfmprd_dir,laps_data_root

        integer max_fcst_times
        parameter (max_fcst_times=200)

	integer ensemble_len,len_dir,i,ct,ctpz
        character EXT*31, directory*255, c_model*30

        character*10  units_2d
        character*125 comment_2d
        character*3 var_3d
        character*3 var_2d
        character*9 a9time,a9time_valid
        character*150 ensemble_dir,enst_dir,ensm_dir
        character*2 domnum_str
!        character*10 mtype

        character*10 ext_fcst_a(n_fields),ext_fcst_a_2d(n_fields)
        character*10 ext_anal_a(n_fields),ext_anal_a_2d(n_fields)
        character*10 var_a(n_fields), var_a_2d(n_fields)
        character*2 c2_region

        character*250 output_dir,cdl_dir,output_file
        character*250 stlaps

!       Specify what is being verified
        data ext_fcst_a /'fua','fua','fua','fua'/ ! 3-D
        data ext_anal_a /'lt1','lw3','lw3','lt1'/
        data var_a      /'ht','u3','v3','t3'/

        data ext_fcst_a_2d /'fsf','fsf','fsf','fsf'/ ! 2-D
        data ext_anal_a_2d /'lsx','lsx','l1s','lsx'/
        data var_a_2d      /'tsf','dsf','r01','fwx'/

        integer n_models_read,lencom,n1,n2
        character*132  comment
        character*3 c_models_read
	character c1, c2



        i4_initial = i4time
        lun_in = 21


	allocate(name3d(n_fields*NZ_L),units3d(n_fields*NZ_L)  
     1   ,lvltype3d(n_fields*NZ_L),com3d(n_fields*NZ_L)
     1   ,lvls3d(n_fields*NZ_L))

	allocate(name2d(n_fields),units2d(n_fields)  
     1   ,lvltype2d(n_fields),com2d(n_fields)
     1   ,lvls2d(n_fields))




!!!!!! test input variables!!!!!!!!!!!!!!!!!!!!!
        n_fcst_times = 38 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!       Get fdda_model_source from static file
        call get_fdda_model_source(c_fdda_mdl_src,n_fdda_models,istatus)

        write(6,*)' n_fdda_models = ',n_fdda_models
        write(6,*)' c_fdda_mdl_src = '
     1            ,(c_fdda_mdl_src(m),m=1,n_fdda_models)

        do itime_fcst = 0,n_fcst_times

          ct=1

          i4_valid = i4_initial + itime_fcst * model_cycle_time_sec

          call make_fnam_lp(i4_valid,a9time_valid,istatus)

          write(6,*)
          write(6,*)
          write(6,*)' Processing time ',itime_fcst,' ',a9time_valid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         3D fields
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          do ifield = 1,n_fields

           write(6,*)

	   mean_fcst_3d=0.

           comment = ' '
           len_com = 1

           var_3d = var_a(ifield)
           call s_len(var_3d,lenvar)


           n_models_read = 0
           do imodel = 1,n_fdda_models

            c_model = c_fdda_mdl_src(imodel)


            if(c_model(1:3) .ne. 'lga')then

              write(6,*)' Processing model ',c_model

              call s_len(c_model,len_model)

              call get_directory('fua',ensemble_dir,len_ensemble)
              write(6,*)' ensemble_dir = ',ensemble_dir
              ensm_dir = ensemble_dir(1:len_ensemble)//'mean/'
              call s_len(ensm_dir,len_ensm)

!                 Read 3d forecast fields

                  ext = ext_fcst_a(ifield)
                  call get_directory(ext,directory,len_dir)
                  DIRECTORY=directory(1:len_dir)//c_model(1:len_model)
     1                                          //'/'

                  call get_lapsdata_3d(i4_initial,i4_valid
     1                          ,NX_L,NY_L,NZ_L       
     1                          ,directory,var_3d
     1                          ,units_2d,comment_2d,var_fcst_3d
     1                          ,istatus)

                  write(*,*)'>',var_3d,'<>',var_fcst_3d(1,1,1),'<'
                  if(istatus .ne. 1)then
                       write(6,*)' Error reading 3D Forecast',c_model
                  else 
                    mean_fcst_3d=mean_fcst_3d+var_fcst_3d


                    write(6,*)'min and max vals for '
     1                                  ,var_a(ifield),c_model,': '
     1                                  ,minval(var_fcst_3d(:,:,:))
     1                                  ,maxval(var_fcst_3d(:,:,:))

		    write(*,*)	

                    n_models_read = n_models_read + 1
                    call s_len(comment,lencom)
                  comment = comment(1:lencom)//c_model(1:len_model)//','
                  endif

          endif ! c_model .ne. lga
         enddo ! model

         write(*,*) 'ct filling mean_fcst_3d_model ',ct
         mean_fcst_3d=mean_fcst_3d/float(n_models_read) 
         mean_fcst_3d_model(:,:,ct:ct+NZ_L-1)=mean_fcst_3d(:,:,:) 
         call s_len(comment,lencom)
         if (n_models_read .lt. 10) then
           c_models_read = '  '//char(n_models_read +48)
         else if (n_models_read .lt. 100) then
           n1 = n_models_read /10
           n2 = mod(n_models_read, 10)
           c1 = char(n1 + 48)
           c2 = char(n2 + 48)
           c_models_read = ' '//c1//c2
         else
           c_models_read = '>99'
         endif
         comment = comment(1:lencom)//'models='//c_models_read
         com3d(ct:ct+NZ_L-1)=comment       
         ct=ct+NZ_L
       
         write(6,*)'min and max vals (mean_fcst_3d) for ' 
     1               ,var_a(ifield),': '
     1               ,minval(mean_fcst_3d(:,:,:))
     1               ,maxval(mean_fcst_3d(:,:,:))

 
        enddo ! ifields


	write(*,*)'n_fdda_models',n_fdda_models

! Write out the 3D stuff using LAPS library routine

        call get_directory('cdl',cdl_dir,len_cdl)
        cdl_dir = cdl_dir(1:len_cdl)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        ct=1
        name3d(ct:ct+NZ_L-1)='HT ' 
        lvls3d(ct:ct+NZ_L-1)=nint(pressures(1:NZ_L)*0.01) 
        j = 1
        do i=NZ_L, 1, -1
          cdl_levels(j) = float(lvls3d(i))
          j = j + 1
        enddo
        ct=ct+NZ_L

!	name3d(ct:ct+NZ_L-1)='RH3'
!        lvls3d(ct:ct+NZ_L-1)=nint(pressures(1:NZ_L)*0.01)
!        ct=ct+NZ_L

!	name3d(ct:ct+NZ_L-1)='SH ' 
!	lvls3d(ct:ct+NZ_L-1)=nint(pressures(1:NZ_L)*0.01)
!	ct=ct+NZ_L

	name3d(ct:ct+NZ_L-1)='U3 '
	lvls3d(ct:ct+NZ_L-1)=nint(pressures(1:NZ_L)*0.01)
	ct=ct+NZ_L

	name3d(ct:ct+NZ_L-1)='V3 '
        lvls3d(ct:ct+NZ_L-1)=nint(pressures(1:NZ_L)*0.01)
        ct=ct+NZ_L

        name3d(ct:ct+NZ_L-1)='T3 '
        lvls3d(ct:ct+NZ_L-1)=nint(pressures(1:NZ_L)*0.01)
        ct=ct+NZ_L


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	write(6,*)'min and max vals - mean_fcst_3d_model: '
     1   ,minval(mean_fcst_3d_model(:,:,:))
     1   ,maxval(mean_fcst_3d_model(:,:,:))
	

       	print*,'Writing LAPS 3D (fua) netcdf file.'
       	print*,'Output directory is: ',ensm_dir(1:len_ensm)
       	print*,'Output directory is: ',trim(ensm_dir)
       	print*,'CDL Directory is: ',trim(cdl_dir)
        call write_laps_lfm(i4_initial,i4_valid,trim(ensm_dir)
     1                ,trim(cdl_dir),'fua'                                          
     1                ,NX_L,NY_L,NZ_L*n_fields,NZ_L*n_fields
     1                ,name3d,lvls3d
     1                ,lvltype3d,units3d,com3d,NZ_L,cdl_levels          
     1                ,mean_fcst_3d_model,istatus)


	  if (istatus /= 1) then
     		print*,'Error writing LAPS 3D (fua) netcdf file.'
 	  else
     		print*,'Done writing 3d netcdf file.'
  	  endif

        enddo ! itime_fcst

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!       2D fields
	write(*,*)'We are starting to deal with 2D fields!!!!!!!!!!!!!'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        call get_directory('fsf',ensemble_dir,len_ensemble)
        write(6,*)' ensemble_dir = ',ensemble_dir
        ensm_dir = ensemble_dir(1:len_ensemble)//'mean/'
        call s_len(ensm_dir,len_ensm)


        do itime_fcst = 0,n_fcst_times

	  ct=1
 
          i4_valid = i4_initial + itime_fcst * model_cycle_time_sec

          call make_fnam_lp(i4_valid,a9time_valid,istatus)

!	  Read 2d forecast fields
          do ifield = 1,n_fields

           comment = ' '
           len_com = 1

           var_2d = var_a_2d(ifield)
           call s_len(var_2d,lenvar)

	   mean_fcst_2d=0.

           n_models_read = 0
           do imodel = 1,n_fdda_models

            c_model = c_fdda_mdl_src(imodel)

            if(c_model(1:3) .ne. 'lga')then

              write(6,*)' Processing model ',c_model

              call s_len(c_model,len_model)


                  ext = ext_fcst_a_2d(ifield)
                  call get_directory(ext,directory,len_dir)
                  DIRECTORY=directory(1:len_dir)//c_model(1:len_model)
     1                                          //'/'


                  call get_lapsdata_2d(i4_initial,i4_valid,directory
     1                  ,var_2d,units_2d,comment_2d,NX_L,NY_L
     1                  ,var_fcst_2d,istatus)

                  write(*,*)'>',var_2d,'<>',var_fcst_2d(1,1),'<'
                  if(istatus .ne. 1)then
                       write(6,*)' Error reading 2D Forecast'
                  else
                    mean_fcst_2d=mean_fcst_2d+var_fcst_2d
                    n_models_read = n_models_read + 1
                    call s_len(comment,lencom)
                  comment = comment(1:lencom)//c_model(1:len_model)//','
                    call s_len(comment,lencom)
                  endif

          endif ! c_model .ne. lga

         enddo ! model

         write(*,*) 'ct filling mean_fcst_2d_model ',ct
         mean_fcst_2d=mean_fcst_2d/float(n_models_read) 
         mean_fcst_2d_model(:,:,ct)=mean_fcst_2d(:,:) 
         call s_len(comment,lencom)
         if (n_models_read .lt. 10) then
           c_models_read = '  '//char(n_models_read +48)
         else if (n_models_read .lt. 100) then
           n1 = n_models_read /10
           n2 = mod(n_models_read, 10)
           c1 = char(n1 + 48)
           c2 = char(n2 + 48)
           c_models_read = ' '//c1//c2
         else
           c_models_read = '>99'
         endif
         comment = comment(1:lencom)//'models='//c_models_read
         com2d(ct:ct)=comment
         ct=ct+1
        
        enddo ! fields

! Write out the 2D stuff using LAPS library routine

        call get_directory('cdl',cdl_dir,len_cdl)
        cdl_dir = cdl_dir(1:len_cdl)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ct=1
        name2d(ct:ct)='TSF'
        lvls2d(ct:ct)=0
        ct=ct+1

        name2d(ct:ct)='DSF'
        lvls2d(ct:ct)=0
        ct=ct+1

        name2d(ct:ct)='R01'
        lvls2d(ct:ct)=0
        ct=ct+1

        name2d(ct:ct)='FWX'
        lvls2d(ct:ct)=0
        ct=ct+1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call write_laps_lfm(i4_initial,i4_valid,trim(ensm_dir)
     1                ,trim(cdl_dir),'fsf'
     1                ,NX_L,NY_L,n_fields,n_fields
     1                ,name2d,lvls2d
     1                ,lvltype2d,units2d,com2d,1,0.
     1                ,mean_fcst_2d_model,istatus)


          if (istatus /= 1) then
                print*,'Error writing LAPS 2D (fsf) netcdf file.'
          else
                print*,'Done writing 2d netcdf file.'
          endif

        enddo ! itime_fcst


 999    write(6,*)' End of subroutine ensemble_mean'




        return

        end


