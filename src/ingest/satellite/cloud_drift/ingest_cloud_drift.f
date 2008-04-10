
!     Steve Albers      Jan-1998        

      integer max_paths
      parameter(max_paths=10)

!     Input file 
      character*170 filename_in
      character*150 dir_in, path_to_cloud_drift(max_paths)
      character*20  cloud_drift_format(max_paths)
      character*255 c_filespec
      integer max_files
      parameter(max_files = 3000)
      character*255 c_fnames(max_files)
      integer i4times(max_files)
 
      character*9 a9_time, a10_to_a9
      character*13 wfo13_time
      character*9 wfo_fname13_to_fname9
      character*10 a10_time

!     Output file
      character*31    ext
      integer       len_dir

      call get_systime(i4time_sys,a9_time,istatus)
      if(istatus .ne. 1)go to 999

!     i4time = (i4time/3600) * 3600

      call get_grid_dim_xy(NX_L,NY_L,istatus)
      if (istatus .ne. 1) then
          write (6,*) 'Error getting horizontal domain dimensions'
          go to 999
      endif
 
      call get_laps_cycle_time(ilaps_cycle_time,istatus)
      if(istatus .eq. 1)then
          write(6,*)' ilaps_cycle_time = ',ilaps_cycle_time
      else
          write(6,*)' Error getting laps_cycle_time'
          stop
      endif

      call get_windob_time_window('CDW',i4_window,istatus)
      if(istatus .ne. 1)stop

      call get_cloud_drift_parms(n_paths_drift
     1                          ,path_to_cloud_drift
     1                          ,cloud_drift_format
     1                          ,istatus)
      if(istatus .ne. 1)stop

      do ipath = 1,n_paths_drift
          dir_in = path_to_cloud_drift(ipath)

          call s_len(dir_in,len_dir_in)
          if(len_dir_in .gt. 0)then
              write(6,*)' path for cloud drift winds = '
     1                   ,dir_in(1:len_dir_in)      
              write(6,*)' data format = ',cloud_drift_format(ipath)
          else
              write(6,*)' Warning: no cloud_drift_path'
              stop
          endif

          if(cloud_drift_format(ipath) .eq. 'CWB_SATOB')then
              c_filespec = dir_in(1:len_dir_in)//'/satob*.dat'
          elseif(cloud_drift_format(ipath) .eq. 'CWB_HDSW')then
              c_filespec = dir_in(1:len_dir_in)//'/hdsw*.dat'
          elseif(cloud_drift_format(ipath) .eq. 'MADIS')then
              c_filespec = dir_in(1:len_dir_in)//'/*'
          else 
              c_filespec = dir_in(1:len_dir_in)
          endif

          call get_file_names(c_filespec,i_nbr_files_ret,c_fnames
     1                       ,max_files,istatus)

!         Open output CDW file 
          lun_out = 11
          ext = 'cdw'
          if(i_nbr_files_ret .gt. 0)then
              write(6,*)' Found ',i_nbr_files_ret,' files in '
     1                 ,c_filespec      
          endif

          if(cloud_drift_format(ipath) .eq. 'NESDIS')then
              i4_contains_early = 0
              i4_contains_late  = 1800

!             Obtain file times from file names
              do i = 1,i_nbr_files_ret
                  call s_len(c_fnames(i),len_fname)
                  call get_directory_length(c_fnames(i),len_dir)
                  a10_time = c_fnames(i)(len_dir+2:len_fname)

                  a9_time = a10_to_a9(a10_time,istatus)
                  if(istatus .ne. 1)then
                      write(6,*)' Bad value for a10_time ',a10_time
                      stop
                  endif

                  call i4time_fname_lp(a9_time,i4times(i),istatus)
                  write(6,*)c_fnames(i)(1:len_fname),i4times(i)
              enddo ! i

          elseif(cloud_drift_format(ipath) .eq. 'AFWA')then
              i4_contains_early = 0
              i4_contains_late = 3600

              call get_file_times(c_filespec,max_files,c_fnames
     1                           ,i4times,i_nbr_files_ret
     1                           ,istatus)

          elseif(cloud_drift_format(ipath) .eq. 'CWB_SATOB')then
              i4_contains_early = 21600
              i4_contains_late = 21600

!             Obtain file times from file names
              do i = 1,i_nbr_files_ret
                  call s_len(c_fnames(i),len_fname)
                  a10_time = c_fnames(i)(len_fname-11:len_fname-4)//'00'       

                  a9_time = a10_to_a9(a10_time,istatus)
                  if(istatus .ne. 1)then
                      write(6,*)' Bad value for a10_time ',a10_time
                      stop
                  endif

                  call i4time_fname_lp(a9_time,i4times(i),istatus)
                  write(6,*)c_fnames(i)(1:len_fname),i4times(i)
              enddo ! i

          elseif(cloud_drift_format(ipath) .eq. 'CWB_HDSW')then
              i4_contains_early = 21600
              i4_contains_late = 21600

!             Obtain file times from file names
              do i = 1,i_nbr_files_ret
                  call s_len(c_fnames(i),len_fname)
                  call get_directory_length(c_fnames(i),len_dir)
                  a10_time = c_fnames(i)(len_fname-11:len_fname-4)//'00'       

                  a9_time = a10_to_a9(a10_time,istatus)
                  if(istatus .ne. 1)then
                      write(6,*)' Bad value for a10_time ',a10_time
                      stop
                  endif

                  call i4time_fname_lp(a9_time,i4times(i),istatus)
                  write(6,*)c_fnames(i)(1:len_fname),i4times(i)
              enddo ! i
           elseif(cloud_drift_format(ipath) .eq. 'MADIS')then
              i4_contains_early = 1800 
              i4_contains_late =  1800
!             Obtain file times from file names
              do i = 1,i_nbr_files_ret
                  call s_len(c_fnames(i),len_fname)
                  wfo13_time=c_fnames(i)(len_fname-12:len_fname)
                  a9_time = wfo_fname13_to_fname9(wfo13_time)
                  call i4time_fname_lp(a9_time,i4times(i),istatus)
                  write(6,*)c_fnames(i)(1:len_fname),i4times(i)
              enddo ! i
          else
              write(6,*)' ERROR, unknown cloud_drift_format '
     1                 ,cloud_drift_format(ipath)

          endif ! cloud_drift_format

!         Loop through ASCII E/W files and choose ones in time window
          write(6,*)' # of data files = ',i_nbr_files_ret
          do i = 1,i_nbr_files_ret
              call make_fnam_lp(i4times(i),a9_time,istatus)
              filename_in = c_fnames(i)

!             Test whether we want the NetCDF file for this time
              i4time_file_earliest = i4time_sys - i4_window
     1                                          - i4_contains_late
              i4time_file_latest =   i4time_sys + i4_window
     1                                          + i4_contains_early
          
              if(i4times(i) .lt. i4time_file_earliest)then
                  write(6,*)' File is too early ',a9_time,i
              elseif(i4times(i) .gt. i4time_file_latest)then
                  write(6,*)' File is too late ',a9_time,i
              else
                  write(6,*)' File is in time window ',a9_time,i
                  write(6,*)' Input file ',filename_in

                  if(cloud_drift_format(ipath) .eq. 'NESDIS')then
!                     Read from the ASCII pirep file

                      call open_ext(lun_out,i4time_sys,ext(1:3),istatus)       
                      if(istatus .ne. 1)then
                          write(6,*)' Error opening output file ',ext
                          stop
                      endif

!                     Write to the opened CDW file
                      call get_cloud_drift_data(i4time_sys,i4_window
     1                                          ,NX_L,NY_L
     1                                          ,filename_in,istatus)

                  elseif(cloud_drift_format(ipath) .eq. 'AFWA')then
                      call open_ext(lun_out,i4time_sys,ext(1:3),istatus)
                      if(istatus .ne. 1)then
                          write(6,*)' Error opening output file ',ext
                          stop
                      endif

                      call get_cloud_drift_afwa(i4time_sys,i4_window
     1                                          ,NX_L,NY_L
     1                                          ,filename_in,istatus)

                  elseif(cloud_drift_format(ipath) .eq. 'CWB_SATOB')then
                      call get_cloud_drift_cwb_satob(i4time_sys
     1                                          ,i4_window,NX_L,NY_L
     1                                          ,filename_in,istatus)

                  elseif(cloud_drift_format(ipath) .eq. 'CWB_HDSW')then
                      call get_cloud_drift_cwb_hdsw(i4time_sys
     1                                          ,i4_window,NX_L,NY_L
     1                                          ,filename_in,istatus)
                  elseif(cloud_drift_format(ipath) .eq. 'MADIS')then
                     call get_cloud_drift_madis(i4time_sys,
     +                  i4_window,filename_in,istatus)
                  else
                      write(6,*)' ERROR, unknown cloud_drift_format ' 
     1                          ,cloud_drift_format(ipath)

                  endif ! cloud_drift_format
              endif ! i4times(i)
          enddo ! i
      enddo ! ipath

 990  close(lun_out) ! Output CDW file

 999  continue

      write(6,*)' End of cloud_drift ingest program'

      end
 
 
       subroutine get_cloud_drift_parms(n_paths_drift
     1                                 ,path_to_cloud_drift
     1                                 ,cloud_drift_format,istatus)

       integer max_paths
       parameter(max_paths=10)

       character*(150) path_to_cloud_drift(max_paths)
       character*(20) cloud_drift_format(max_paths)

       namelist /cloud_drift_nl/ n_paths_drift, path_to_cloud_drift
     1                                        , cloud_drift_format
 
       character*150 static_dir,filename
 
       call get_directory('nest7grid',static_dir,len_dir)

       filename = static_dir(1:len_dir)//'/cloud_drift.nl'
 
       open(1,file=filename,status='old',err=900)
       read(1,cloud_drift_nl,err=901)
       close(1)

       istatus = 1
       return

  900  print*,'error opening file ',filename
       istatus = 0
       return

  901  print*,'error reading cloud_drift_nl in ',filename
       write(*,cloud_drift_nl)
       istatus = 0
       return
 
       end

