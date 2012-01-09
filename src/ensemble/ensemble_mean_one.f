      program HMT_ensemble_mean_main

! Program calculates mean for a single time

! Usage: ensemble_mean.exe [h]hh LAPS_DATA_ROOT

! Program calculates ensmble mean as well as time mean
! Isidora Jankov, Linda Wharton and Steve Albers July 2010
! modified for HMT 2011-12 LW Nov 2011

      character*9 a9time
      character*80 buffer
      character*300 laps_data_root,dir_t,filenamet
      character*250 stlaps
      integer      NX_L,NY_L,NZ_L,i4time,i,istatus
      integer      model_cycle_time_sec, len_dir_t
      integer      run_hr, debug
      real         r_missing_data
      integer, parameter :: nprmax=150
      real, dimension(nprmax) :: pressures
      integer, parameter :: lun=120

       include '../include/lapsparms.for'

! Begin

      debug = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get parameters from laps_data_root
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Parameter for model cycle interval
      model_cycle_time_sec = 10800

      call getarg(1,buffer)
      read(buffer,*) run_hr
      call getarg(2,laps_data_root)
      print*, 'run_hr >',run_hr,'<'
      print*, 'l_d_r >',trim(laps_data_root),'<'

!      call get_systime(i4time,a9time,istatus)
!      if(istatus .ne. 1)go to 999
!      write(6,*)' systime = ',a9time

! this is written by something in /lfs0/projects/hmtb/dwr_domains/laps/d01/time
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

      print*, 'LW i4time, a9time ',i4time, a9time

      call get_grid_dim_xy(NX_L,NY_L,istatus)
      if (istatus .ne. 1) then
        write (6,*) 'Error getting horizontal domain dimensions'
        go to 999
      else
        if (debug .eq. 1) then
          print*,'NX_L = ',NX_L
          print*,'NY_L = ',NY_L
        endif
      endif

      call get_laps_dimensions(NZ_L,istatus)
      if (istatus .ne. 1) then
        write (6,*) 'Error getting vertical domain dimension'
        go to 999
      else
        if (debug .eq. 1) then
          print*,'NZ_L = ',NZ_L
        endif
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
     1                   NX_L,NY_L,
     1                   NZ_L,
     1                   pressures,
     1                   run_hr,
     1                   r_missing_data,
     1                   istatus)

999   continue

      end

!******************************************************************

      subroutine ensemble_mean(i4time,model_cycle_time_sec,
     1                  NX_L,NY_L,
     1                  NZ_L,
     1                  pressures,
     1                  run_hr,
     1                  r_missing_data,
     1                  istatus)

       include '../include/lapsparms.for'

      integer      i4time,istatus,debug
      integer      model_cycle_time_sec, run_hr
      integer,intent(in) ::   NX_L,NY_L,NZ_L
      real         pressures(NZ_L)
      real         r_missing_data

      character(len=15),  allocatable, dimension(:) :: ensDir
      character(len=31),  allocatable, dimension(:) :: ext
      character(len=3),   allocatable, dimension(:) :: var
      integer,            allocatable, dimension(:) :: lvl2d

      integer n_2d_fields, calcmaxes, n_2d_fields_write
      real, allocatable, dimension(:,:,:) :: mean_fcst_2d_model !(NX_L,NY_L,n_2d_fields)
      real, allocatable, dimension(:,:) :: mean_fcst_2d       !(NX_L,NY_L)
      real, allocatable, dimension(:,:) :: var_fcst_2d        !(NX_L,NY_L)
      real, allocatable, dimension(:,:) :: r01_max_2d         !(NX_L,NY_L)
      real, allocatable, dimension(:,:) :: rto_max_2d         !(NX_L,NY_L)
      character(len=3),   allocatable, dimension(:) :: name2d    !n_2d_fields
      character(len=10),  allocatable, dimension(:) :: units2d   !n_2d_fields
      character(len=4),   allocatable, dimension(:) :: lvltype2d !n_2d_fields
      character(len=132), allocatable, dimension(:) :: com2d     !n_2d_fields
      integer, allocatable, dimension(:) :: lvls2d               !n_2d_fields

      integer n_3d_fields
      real, allocatable, dimension(:,:,:) :: mean_fcst_3d_model !(NX_L,NY_L,NZ_L*n_3d_fields)
      real, allocatable, dimension(:,:,:) :: mean_fcst_3d       !(NX_L,NY_L,NZ_L)
      real, allocatable, dimension(:,:,:) :: var_fcst_3d        !(NX_L,NY_L,NZ_L)
      character(len=3),   allocatable, dimension(:) :: name3d
      character(len=10),  allocatable, dimension(:) :: units3d
      character(len=4),   allocatable, dimension(:) :: lvltype3d
      character(len=132), allocatable, dimension(:) :: com3d
      integer, allocatable, dimension(:) :: lvls3d

      real cdl_levels(NZ_L)

!      logical write_to_lapsdir

      integer       i4_initial, i4_valid
!      integer       maxbgmodels
!      parameter     (maxbgmodels=10)
      character*30  c_fdda_mdl_src(maxbgmodels)
      character*300 lfmprd_dir,laps_data_root

      integer max_fcst_times
      parameter (max_fcst_times=200)

      integer ensemble_len,len_dir,ct,ctpz
      character directory*255, c_model*30

      character*3   var_2d
      character*10  units_2d
      character*125 comment_2d
      character*4   lvltype_2d
      character*3   var_3d
      character*10  units_3d
      character*125 comment_3d
      character*4   lvltype_3d
      character*9 a9time, a9time_valid
      character*150 ensemble_dir, mean_dir, base_directory
      character*2 domnum_str
!      character*10 mtype

      character*2 c2_region

      character*250 output_dir,cdl_dir,output_file
      character*250 stlaps
      character*120 file_name
      character*9 gtime

      integer n_models_read,lencom,n1,n2,len,numEnsMemb,i,j
      integer nmr_2d, nmr_3d
      character*132  comment
      character*3 c_models_read
      character c1, c2
      character*5       fcst_hh_mm
      integer     ext_len, fn_length
      integer     no_nmm_umf, no_nmm_s8a

! Begin
      debug = 0

      print*, 'LW inside run_hr >',run_hr,'< i4time ',i4time

! members to skip for now
! Skip nmm models for mean of umf
      no_nmm_umf = 1
! Skip nmm models for mean of s8a
      no_nmm_s8a = 0


! LW Stubs put in for driving by namelist
!   namelist format:
!   base directory
!   number ens members
!   ensDir(1)
!   :
!   ensDir(numEnsMemb)
!   n_fcst_times
!  !ext  var  lvl
!   fua  u3  925
!   fua  v3  925
!   fsf  tsf 0
!   fsf  dsf 0
!!!!!! Hardwired for now LW !!!!!!!!!!!!!!!!!!!
      mean_dir =
     1'/lfs0/projects/hmtb/dwr_domains/laps/d01/lapsprd/fsf/mean/'
      print*, '==================================================='
      print*, 'Hardwired mean_dir=',trim(mean_dir)
      numEnsMemb = 9
      allocate(ensDir(numEnsMemb))
      ensDir(1) = 'arw-tom-gep0'
      ensDir(2) = 'arw-fer-gep1'
      ensDir(3) = 'arw-sch-gep2'
      ensDir(4) = 'arw-tom-gep3'
      ensDir(5) = 'nmm-fer-gep4'
      ensDir(6) = 'arw-fer-gep5'
      ensDir(7) = 'arw-sch-gep6'
      ensDir(8) = 'arw-tom-gep7'
      ensDir(9) = 'nmm-fer-gep8'
!      n_fcst_times = 38
      n_fcst_times = 1 !process ONLY the hour passed in
      n_3d_fields = 0
      n_2d_fields = 17
      allocate(ext(n_2d_fields))
      allocate(var(n_2d_fields))
      allocate(lvl2d(n_2d_fields))
      ext(1) = 'fua'
      var(1) = 'u3 '
      lvl2d(1) = 925
      ext(2) = 'fua'
      var(2) = 'v3 '
      lvl2d(2) = 925
      ext(3) = 'fsf'
      var(3) = 'tsf'
      lvl2d(3) = 0
      ext(4) = 'fsf'
      var(4) = 'dsf'
      lvl2d(4) = 0
      ext(5) = 'fsf'
      var(5) = 'rh '
      lvl2d(5) = 0
      ext(6) = 'fsf'
      var(6) = 'usf'
      lvl2d(6) = 0
      ext(7) = 'fsf'
      var(7) = 'vsf'
      lvl2d(7) = 0
      ext(8) = 'fsf'
      var(8) = 'r01'
      lvl2d(8) = 0
      ext(9) = 'fsf'
      var(9) = 'rto'
      lvl2d(9) = 0
      ext(10) = 'fsf'
      var(10) = 's01'
      lvl2d(10) = 0
      ext(11) = 'fsf'
      var(11) = 'sto'
      lvl2d(11) = 0
      ext(12) = 'fsf'
      var(12) = 'umf'
      lvl2d(12) = 0
      ext(13) = 'fsf'
      var(13) = 'tw0'
      lvl2d(13) = 0
      ext(14) = 'fsf'
      var(14) = 'tw1'
      lvl2d(14) = 0
      ext(15) = 'fsf'
      var(15) = 's8a'
      lvl2d(15) = 0
      ext(16) = 'fsf'
      var(16) = 'vis'
      lvl2d(16) = 0
      ext(17) = 'fsf'
      var(17) = 'fwx'
      lvl2d(17) = 0

! LW generate grid that has the max value at each gridpoint for r01 and rto grids
      calcmaxes = 1
      print*, '==================================================='
      print*, 'Option to output fields containing Max Incremental '//
     1'Precip'
      print*, 'and Max Model Run Total Precip Accum at each '//
     1'gridpoint'
      if (calcmaxes .eq. 1) then
        print*, 'is turned ON'
      else
        print*, 'is turned OFF'
      endif
      print*, '==================================================='
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Allocate variables to hold data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (n_3d_fields .gt. 0) then ! Allocate 3d variables
        allocate(name3d(n_3d_fields*NZ_L))
        allocate(units3d(n_3d_fields*NZ_L))
        allocate(lvltype3d(n_3d_fields*NZ_L))
        allocate(com3d(n_3d_fields*NZ_L))
        allocate(lvls3d(n_3d_fields*NZ_L))
        allocate(mean_fcst_3d_model(NX_L,NY_L,NZ_L*n_3d_fields))
        allocate(mean_fcst_3d(NX_L,NY_L,NZ_L))
        allocate(var_fcst_3d(NX_L,NY_L,NZ_L))
      endif

      if (n_2d_fields .gt. 0) then ! Allocate 2d variables
        if (calcmaxes .eq. 1) then
          n_2d_fields_write = n_2d_fields + 2
          allocate(r01_max_2d(NX_L,NY_L))
          allocate(rto_max_2d(NX_L,NY_L))
        else
          n_2d_fields_write = n_2d_fields
        endif
        allocate(mean_fcst_2d(NX_L,NY_L))
        allocate(var_fcst_2d(NX_L,NY_L))
        allocate(name2d(n_2d_fields_write))
        allocate(units2d(n_2d_fields_write))
        allocate(lvltype2d(n_2d_fields_write))
        allocate(com2d(n_2d_fields_write))
        allocate(lvls2d(n_2d_fields_write))
        allocate(mean_fcst_2d_model(NX_L,NY_L,n_2d_fields_write))
      endif

! LW old way...get model extensions from nest7grid.parms variable FDDA_MODEL_SOURCE
!      Get fdda_model_source from static file
!      call get_fdda_model_source(c_fdda_mdl_src,n_fdda_models,istatus)

      write(6,*)' Number Ensemble Members = ',numEnsMemb
      print*, '---------------------------------------------------'
      write(6,*)' Ensembles = '
     1          ,(ensDir(m),m=1,numEnsMemb)
      print*, '---------------------------------------------------'
      write(6,*)' Directories = '
      do i = 1, numEnsMemb
        call s_len(ensDir(i),len_model)
        call get_directory('fsf',base_directory,len_dir)
        directory=base_directory(1:len_dir)//ensDir(i)(1:len_model)
        write(6,*)'   ',trim(directory)
      enddo

      call get_directory('cdl',cdl_dir,len_cdl)
      print*, '---------------------------------------------------'
      print*, 'cdl_dir ',trim(cdl_dir)
      call s_len(mean_dir,len_ensm)
      print*, '---------------------------------------------------'
      print*, 'Output mean_dir  ',trim(mean_dir)

      i4_initial = i4time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Loop over fcst times, writing out mean file for each time
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      itime_fcst = run_hr

        ct=1

!        i4_valid = i4_initial + itime_fcst * model_cycle_time_sec
        i4_valid = i4_initial + run_hr * 3600
        print*, i4_initial, run_hr

        call make_fnam_lp(i4_valid,a9time_valid,istatus)

        write(6,*)'=================================================='
100     format(a21,i4,a12,a9)
        write(6,100)' Processing: fcst hr=',run_hr,
     1'  fcst time=',a9time_valid
        write(6,*)'=================================================='
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Process 3d fields
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (n_3d_fields .gt. 0) then !LW old logic

        do ifield = 1,n_3d_fields

          write(*,*)'========= Starting to process 3d fields =========='

          mean_fcst_3d=0.

          comment = ' '
          len_com = 1

          var_3d = var(ifield)
          call s_len(var_3d,lenvar)

          n_models_read = 0
          do imodel = 1,n_fdda_models

            c_model = c_fdda_mdl_src(imodel)


            if(c_model(1:3) .ne. 'lga')then

              write(6,*)' Processing model ',c_model

              call s_len(c_model,len_model)

              call get_directory('fua',ensemble_dir,len_ensemble)
              write(6,*)' ensemble_dir = ',ensemble_dir
              mean_dir = ensemble_dir(1:len_ensemble)//'mean/'
              call s_len(mean_dir,len_ensm)

!             Read 3d forecast fields

              call get_directory(ext(ifield),directory,len_dir)
              DIRECTORY=directory(1:len_dir)//c_model(1:len_model)
     1                                          //'/'

              call get_lapsdata_3d(i4_initial,i4_valid
     1                          ,NX_L,NY_L,NZ_L
     1                          ,directory,var_3d
     1                          ,units_2d,comment_2d,var_fcst_3d
     1                          ,istatus)

              write(*,*)'>',var_3d,'<>',var_fcst_3d(1,1,1),'<'
              if(istatus .ne. 1)then
                write(6,*)' Error reading 3d Forecast',c_model
              else
                mean_fcst_3d=mean_fcst_3d+var_fcst_3d

                write(6,*)'min and max vals for '
!     1                   ,var_a(ifield),c_model,': '
     1                   ,minval(var_fcst_3d(:,:,:))
     1                   ,maxval(var_fcst_3d(:,:,:))

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
!     1               ,var_a(ifield),': '
     1               ,minval(mean_fcst_3d(:,:,:))
     1               ,maxval(mean_fcst_3d(:,:,:))


        enddo ! ifields


        write(*,*)'n_fdda_models',n_fdda_models

! Write out the 3d stuff using LAPS library routine

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

!       name3d(ct:ct+NZ_L-1)='RH3'
!        lvls3d(ct:ct+NZ_L-1)=nint(pressures(1:NZ_L)*0.01)
!        ct=ct+NZ_L

!       name3d(ct:ct+NZ_L-1)='SH '
!       lvls3d(ct:ct+NZ_L-1)=nint(pressures(1:NZ_L)*0.01)
!       ct=ct+NZ_L

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


        print*,'Writing LAPS 3d (fua) netcdf file.'
        print*,'Output directory is: ',mean_dir(1:len_ensm)
        print*,'Output directory is: ',trim(mean_dir)
        print*,'CDL Directory is: ',trim(cdl_dir)
        call write_laps_lfm(i4_initial,i4_valid,trim(mean_dir)
     1                ,trim(cdl_dir),'fua'
     1                ,NX_L,NY_L,NZ_L*n_fields,NZ_L*n_fields
     1                ,name3d,lvls3d
     1                ,lvltype3d,units3d,com3d,NZ_L,cdl_levels
     1                ,mean_fcst_3d_model,istatus)


          if (istatus /= 1) then
                print*,'Error writing LAPS 3d (fua) netcdf file.'
          else
                print*,'Done writing 3d netcdf file.'
          endif

        endif ! if (n_3d_fields .gt. 0)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       2d fields
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (n_2d_fields .gt. 0) then

          write(*,*)'========= Starting to process 2d fields =========='

          ! reset vars to zero
          mean_fcst_2d_model = 0.0
          name2d = '   '
          units2d = '          '
          lvltype2d = '    '
          com2d = ''
          lvls2d = 0
          r01_max_2d = 0.0
          rto_max_2d = 0.0
          nmr_2d = 0

          do ifield = 1,n_2d_fields  !loop through 2d fields

            ! reset vars to zero
            mean_fcst_2d = 0.0
            var_fcst_2d = 0.0

            if (debug .eq. 1) then
              print*, 'Field ',ifield,' var ',var(ifield)
            endif

            call get_directory(ext(ifield),base_directory,len_dir)

            var_2d = var(ifield)
            call s_len(var_2d,lenvar)
            lvl_2d = lvl2d(ifield)

      print*, '==================================================='
            print*, 'Retrieving data for variable ',var_2d
      print*, '==================================================='

            n_models_read = 0

            do member = 1, numEnsMemb  !loop through ensemble members
              c_model = ensDir(member)
              call s_len(c_model,len_model)
        directory=base_directory(1:len_dir)//c_model(1:len_model)//'/'

              if (debug .eq. 1) then
                call make_fnam_lp(i4_initial,gtime,istatus)
           call make_fcst_time(i4_valid,i4_initial,fcst_hh_mm,istatus)
                call cvt_fname_v3(directory,gtime,fcst_hh_mm,
     1                       ext(ifield),ext_len,
     1                       file_name,fn_length,istatus)

                print*, 'file_name ',trim(file_name)
                print*, 'NX_L ',NX_L, ' NY_L ',NY_L
              endif


! LW temporary construct to exclude nmm models from umf and s8a mean calculations
             if ((no_nmm_umf .eq. 1) .and. (var_2d .eq. 'umf')
     1           .and. ((trim(ensDir(member)) .eq. 'nmm-fer-gep4')
     1           .or. (trim(ensDir(member)) .eq. 'nmm-fer-gep8'))) then
               print*, 'umf for ',trim(ensDir(member)),
     1' EXCLUDED from mean'
             else if ((no_nmm_s8a .eq. 1) .and. (var_2d .eq. 's8a')
     1           .and. ((trim(ensDir(member)) .eq. 'nmm-fer-gep4')
     1           .or. (trim(ensDir(member)) .eq. 'nmm-fer-gep8'))) then
               print*, 'sa8 for ',trim(ensDir(member)),
     1' EXCLUDED from mean'
             else
              call read_laps(i4_initial,i4_valid,directory,ext(ifield),
     1              NX_L,NY_L,1,1,var_2d,lvl_2d,lvltype_2d,
     1              units_2d,comment_2d,var_fcst_2d,istatus)

              if(istatus .ne. 1)then
                write(6,*)' Error reading 2d Forecast ',
     1                     var_2d, c_model
              else
                n_models_read = n_models_read + 1
                if (debug .eq. 1) then
                  write(*,*)'>',var_2d,'<>',var_fcst_2d(1,1),'<'
                endif

                if (calcmaxes .eq. 1) then
! calculate the max of all members at each gridpoint for r01 and rto
                  if (var(ifield) .eq. 'r01') then
                    do i = 1, NX_L
                      do j = 1, NY_L
                         if (var_fcst_2d(i,j) .gt. r01_max_2d(i,j)) then
                           r01_max_2d(i,j) = var_fcst_2d(i,j)
                         endif
                       enddo
                    enddo
                  endif
                  if (var(ifield) .eq. 'rto') then
                    do i = 1, NX_L
                      do j = 1, NY_L
                         if (var_fcst_2d(i,j) .gt. rto_max_2d(i,j)) then
                           rto_max_2d(i,j) = var_fcst_2d(i,j)
                         endif
                       enddo
                    enddo
                  endif
                endif  !if calcmaxes .eq. 1

                write(6,*)'min and max for ',trim(c_model)
     1                    ,minval(var_fcst_2d(:,:))
     1                    ,maxval(var_fcst_2d(:,:))
                mean_fcst_2d=mean_fcst_2d+var_fcst_2d

! the construct for u3 and v3 at 925 mb is to allow them to be written
! out to the 2d fsf file
                if ((var(ifield) .eq. 'u3 ') .and.
     1              (lvl2d(ifield) .eq. 925)) then
                  name2d(ifield) = 'u01'
                  lvls2d(ifield) = 0
                else
                  if ((var(ifield) .eq. 'v3 ') .and.
     1                (lvl2d(ifield) .eq. 925)) then
                    name2d(ifield) = 'v01'
                    lvls2d(ifield) = 0
                  else
                    name2d(ifield) = var(ifield)
                    lvls2d(ifield) = lvl2d(ifield)
                  endif
                endif
                call s_len(units_2d,len)
                units2d(ifield) = units_2d(1:len)
                call s_len(lvltype_2d,len)
                lvltype2d(ifield) = lvltype_2d(1:len)
                call s_len(comment_2d,lencom)
                call s_len(com2d(ifield),lencom)
                com2d(ifield) = com2d(ifield)(1:lencom)
     1//c_model(1:len_model)//','
              endif ! read fcst correctly

              nmr_2d = max(n_models_read, nmr_2d)

             endif
            enddo ! members

            if (debug .eq. 1) then
              print*, 'LW after members >',name2d(ifield),'<>',
     1lvl2d(ifield),'<>',com2d(ifield),'<'
            endif

            if (n_models_read .eq. 0) then !no models to average
              print*, 'No models read for var ',var(ifield)
              print*, 'Setting data values to ',r_missing_data
              mean_fcst_2d = r_missing_data
            else
              if ((no_nmm_umf .eq. 1) .or. (no_nmm_s8a .eq. 1) .or.
     1            (n_models_read .ne. numEnsMemb)) then
                print*, 'Included in mean: n_models_read = ',
     1                   n_models_read
              endif

              write(*,*) 'Filling mean_fcst_2d_model ', var(ifield)
              mean_fcst_2d = mean_fcst_2d / float(n_models_read)
              write(6,*)'min and max of mean for '
     1                  ,var(ifield),':'
     1                  ,minval(mean_fcst_2d(:,:))
     1                  ,maxval(mean_fcst_2d(:,:))
            endif
            mean_fcst_2d_model(:,:,ifield) = mean_fcst_2d(:,:)
            if (n_models_read .lt. 10) then
              c_models_read = '  '//char(n_models_read +48)
            else
              if (n_models_read .lt. 100) then
                n1 = n_models_read /10
                n2 = mod(n_models_read, 10)
                c1 = char(n1 + 48)
                c2 = char(n2 + 48)
                c_models_read = ' '//c1//c2
              else
                c_models_read = '>99'
              endif
            endif
            call s_len(com2d(ifield),lencom)
! LW ensemble names currently too long to add models=c_models_read at end of string
!            com2d(ifield) = com2d(ifield)(1:lencom)//' models='//c_models_read
            if (debug .eq. 1) then
              print*, 'LW c_models_read ',c_models_read
              print*, 'LW com2d(ifield) / lencom ',com2d(ifield) ,lencom
            endif

          enddo ! fields

! Write out the 2d stuff using LAPS library routine

          if (debug .eq. 1) then
            print*,'LW after fields'
          endif

          if (nmr_2d .eq. 0) then
      print*, '==================================================='
            print*, 'No models read...no ensemble mean file output'//
     1' for ', a9time_valid,' time'
      print*, '==================================================='
          else
            call get_directory('cdl',cdl_dir,len_cdl)
            cdl_dir = cdl_dir(1:len_cdl)

            call s_len(mean_dir,len_ensm)

            if (debug .eq. 1) then
              print*, 'LW before write_laps_lfm '
              print*, 'cdl_dir = ',trim(cdl_dir)
              print*,'name2d ',name2d
              print*,'lvls2d ',lvls2d
              print*,'lvltype2d ',lvltype2d
              print*,'units2d ',units2d
              print*,'com2d ',com2d
            endif

            if (calcmaxes .eq. 1) then !add mr0 and mrt to output
              name2d(n_2d_fields+1) = 'mr0'
              lvls2d(n_2d_fields+1) = 0
              lvltype2d(n_2d_fields+1) = 'AGL '
              units2d(n_2d_fields+1) = 'M'
              com2d(n_2d_fields+1) = 'Max incremental precip accum on'//
     1' a point by point basis across ensemble members'
              mean_fcst_2d_model(:,:,n_2d_fields+1) = r01_max_2d(:,:)
      print*, '==================================================='
              write(6,*)'Max Incremental Precip:           '
     1                  ,maxval(r01_max_2d(:,:))
              name2d(n_2d_fields+2) = 'mrt'
              lvls2d(n_2d_fields+2) = 0
              lvltype2d(n_2d_fields+2) = 'AGL '
              units2d(n_2d_fields+2) = 'M'
              com2d(n_2d_fields+2) = 'Max run total precip accum on'//
     1' a point by point basis across ensemble members'
              mean_fcst_2d_model(:,:,n_2d_fields+2) = rto_max_2d(:,:)
      print*, '==================================================='
              write(6,*)'Max Model Run Total Precip Accum: '
     1                  ,maxval(rto_max_2d(:,:))
      print*, '==================================================='
            endif

            call write_laps_lfm(i4_initial,i4_valid,trim(mean_dir)
     1                ,trim(cdl_dir),'fsf',NX_L,NY_L
     1                ,n_2d_fields_write,n_2d_fields_write
     1                ,name2d,lvls2d
     1                ,lvltype2d,units2d,com2d,1,0.
     1                ,mean_fcst_2d_model,istatus)


            if (debug .eq. 1) then
              print*, 'LW after write_laps_lfm / istatus ',istatus
            endif
            if (istatus /= 1) then
              print*,'Error writing ensemble mean 2d (fsf) netcdf file.'
            else
              print*,'Done writing ensemble mean 2d netcdf file.'
            endif
          endif ! if nmr_2d .eq. 0
        endif ! if (n_2d_fields .gt. 0)

      print*, '==================================================='
 999  write(6,*)' End of subroutine ensemble_mean'
      print*, '==================================================='

      return

      end
