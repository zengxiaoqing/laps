
!     Steve Albers      Jan-1998        

      integer max_paths
      parameter(max_paths=2)

!     Input file 
      character*70 filename_in
      character*150 dir_in, path_to_cloud_drift(max_paths)
      character*255 c_filespec
      integer max_files
      parameter(max_files = 3000)
      character*255 c_fnames(max_files)
      integer i4times(max_files)

      character*9 a9_time, a10_to_a9
      character*10 a10_time

!     Output file
      character*31    ext
      integer*4       len_dir

      character*40 c_vars_req
      character*100 c_values_req

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

      i4_window = 7200
      lag_time_report = 1800

!     Open output PIN file to append
      ext = 'cdw'
      call open_lapsprd_file(11,i4time_sys,ext(1:3),istatus)

      call get_cloud_drift_parms(path_to_cloud_drift,istatus)

      do ipath = 1,max_paths
          dir_in = path_to_cloud_drift(ipath)

          call s_len(dir_in,len_dir_in)
          if(istatus .eq. 1)then
              write(6,*)' path for cloud drift winds = '
     1                   ,dir_in(1:len_dir_in)      
          else
              write(6,*)' Error getting cloud_drift_path'
              stop
          endif

          c_filespec = dir_in(1:len_dir_in)

          call get_file_names(c_filespec,i_nbr_files_ret,c_fnames
     1                       ,max_files,istatus)

!         Obtain file times from file names
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

!         Loop through ASCII E/W files and choose ones in time window
          write(6,*)' # of ascii data files = ',i_nbr_files_ret
          do i = 1,i_nbr_files_ret
              call make_fnam_lp(i4times(i),a9_time,istatus)
              filename_in = c_fnames(i)

!             Test whether we want the NetCDF file for this time
              i4time_file_earliest = i4time_sys - i4_window
     1                                          - lag_time_report
              i4time_file_latest =   i4time_sys + i4_window
          
              if(i4times(i) .lt. i4time_file_earliest)then
                  write(6,*)' File is too early ',a9_time,i
              elseif(i4times(i) .gt. i4time_file_latest)then
                  write(6,*)' File is too late ',a9_time,i
              else
                  write(6,*)' File is in time window ',a9_time,i
                  write(6,*)' Input file ',filename_in

!                 Read from the ASCII pirep file and write to the opened PIN file
                  call get_cloud_drift_data(i4time_sys,i4_window
     1                                      ,NX_L,NY_L
     1                                      ,filename_in,istatus)
              endif ! i4times(i)
          enddo ! i
      enddo ! ipath

 990  close(11) ! Output PIN file

 999  continue

      write(6,*)' End of cloud_drift ingest program'

      end
 
 
       subroutine get_cloud_drift_parms(path_to_cloud_drift,istatus)

       integer max_paths
       parameter(max_paths=2)

       character*150 path_to_cloud_drift(max_paths),filename
       namelist /cloud_drift_nl/ path_to_cloud_drift
 
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


        function a10_to_a9(a10_time,istatus)

!       Convert a10_time (yyMMddhhmm) to a9_time (yydddhhmm)

        character*10 a10_time
        character*9 a10_to_a9, a8_to_a9, a9_time
        character*8 a8_time

        a8_time = a10_time(1:8)
        a9_time = a8_to_a9(a8_time) 
        a9_time = a9_time(1:7)//a10_time(9:10)

        a10_to_a9 = a9_time

        istatus = 1

        return
        end


        function a8_to_a9(a8_time)

!       Convert a8_time (yyMMddhh) to a9_time (yydddhhmm)

        character*9 a8_to_a9
        character*9 a9
        character*8 a8_time

        integer*4 imon_a(12)

        data imon_a/0,31,59,90,120,151,181,212,243,273,304,334/

        read(a8_time,1)iyr,imn,idy,ih
1       format(i2,i2,i2,i2)

        id = imon_a(imn)

        idays = id + idy

!       Decide whether to add a day for leap year.
        if(iyr .eq. (iyr / 4) * 4 )then
            if(imn .ge. 3)then
                idays = idays + 1
            endif
        endif

        write(a9,2)iyr,idays,ih,im
2       format(i2,i3,i2,i2)

        if(a9(3:3) .eq. ' ')a9(3:3) = '0'
        if(a9(4:4) .eq. ' ')a9(4:4) = '0'
        if(a9(6:6) .eq. ' ')a9(6:6) = '0'
        if(a9(8:8) .eq. ' ')a9(8:8) = '0'

        a8_to_a9 = a9

        return
        end

