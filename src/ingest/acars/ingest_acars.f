
      subroutine ingest_acars(istatus)

!     Input file 
      character*170 filename_in
      character*9 a9_time
      character*180 dir_in
      character*3 ext_in,fname_fmt
      character*8 c8_file_fmt,c8_project
      character*255 c_filespec
      integer max_files
      parameter(max_files = 3000)
      character*255 c_fnames(max_files)
      integer i4times(max_files)
      character*13 filename13, cvt_i4time_wfo_fname13
      logical l_use_tamdar ! applies to non-WFO, non_RSA runs

!     Output file
      character*31    ext

      character*40 c_vars_req
      character*180 c_values_req

      call get_systime(i4time_sys,a9_time,istatus)
      if(istatus .ne. 1)go to 999

!     i4time = (i4time/3600) * 3600

      call get_grid_dim_xy(NX_L,NY_L,istatus)
      if (istatus .ne. 1) then
          write (6,*) 'Error getting horizontal domain dimensions'
          go to 999
      endif

      call get_c8_project(c8_project,istatus)
      if (istatus .ne. 1) then
          write (6,*) 'Error getting c8_project'
          go to 999
      endif

      call get_l_use_tamdar(l_use_tamdar,istatus)
      if (istatus .ne. 1) then
          write (6,*) 'Error getting l_use_tamdar'
          go to 999
      endif

      lag_time_report = 3600

!     Get path to input files (/public NetCDF format)
      c_vars_req = 'path_to_qc_acars'
      call get_static_info(c_vars_req,c_values_req,1,istatus)
      if(istatus .eq. 1)then
          write(6,*)c_vars_req(1:30),' = ',c_values_req
          dir_in = c_values_req
      else
          write(6,*)' Error getting ',c_vars_req
          return
      endif

      if(c8_project .eq. 'NIMBUS' .or. c8_project .eq. 'CWB')then
          ext_in = 'cdf'
          call s_len(dir_in,len_dir_in)
          c_filespec = dir_in(1:len_dir_in)//'/'//'*00q.'//ext_in

!         Wait for latest input data (only for NIMBUS format)
          i4time_now = i4time_now_gg()
          i4_hour = (i4time_now/3600) * 3600        ! i4time at the top of the hour
          minutes_now = (i4time_now - i4_hour) / 60

          if(minutes_now .ge. 19 .and. minutes_now .lt. 22)then
              i4time_desired = i4_hour
              i4_check_interval = 10
              i4time_stop_waiting = i4_hour + 22*60
              i4_total_wait = min(i4time_stop_waiting - i4time_now, 120)
              i4_thresh_age = 3600

              call wait_for_data(c_filespec,i4time_desired
     1                   ,i4_check_interval,i4_total_wait
     1                   ,i4_thresh_age       ! Only loop through the waiting
                                              ! if data is younger than this
                                              ! threshold
     1                   ,istatus)
          endif ! within time range to wait for data

          c8_file_fmt = 'NIMBUS'

      elseif(c8_project .eq. 'AFWA')then
          ext_in = 'ac'
          c8_file_fmt = 'AFWA'
          call s_len(dir_in,len_dir_in)
          c_filespec = dir_in(1:len_dir_in)//'/'//'*00q.'//ext_in

      else ! Create filespec for WFO format
          ext_in = 'wfo'
          c8_file_fmt = 'WFO'
          call s_len(dir_in,len_dir_in)
          c_filespec = dir_in(1:len_dir_in)

      endif

      fname_fmt = ext_in

!     Get list of files
      call get_file_times(c_filespec,max_files,c_fnames
     1                      ,i4times,i_nbr_files_ret,istatus)

!     Check to see if files are named (presumably) with the WFO convention 
      if(i_nbr_files_ret .eq. 0 .and. fname_fmt .ne. 'wfo')then
          write(6,*)' No raw files with filename convention of '
     1              ,fname_fmt
          write(6,*)
     1        ' Try for wfo filename convention (e.g. MADIS ACARS)...'       
          ext_in = 'wfo'
          fname_fmt = 'wfo'
          call s_len(dir_in,len_dir_in)
          c_filespec = dir_in(1:len_dir_in)
          call get_file_times(c_filespec,max_files,c_fnames
     1                          ,i4times,i_nbr_files_ret,istatus)
      endif

!     Open output PIN file for appending
      if(i_nbr_files_ret .gt. 0 .and. istatus .eq. 1)then
          ext = 'pin'
      else
          write(6,*)' No raw data files identified of fname fmt '
     1              ,fname_fmt
          goto999
      endif

!     Get ACARS Time Window
      call get_windob_time_window('ACARS',i4_wind_ob,istatus)
      if(istatus .ne. 1)goto 997

      call get_tempob_time_window('ACARS',i4_temp_ob,istatus)
      if(istatus .ne. 1)goto 997

      i4_acars_window = max(i4_wind_ob,i4_temp_ob)

!     Loop through ACARS files and choose ones in time window
      write(6,*)' # of raw ACARS files = ',i_nbr_files_ret
      write(6,*)' File Format = ',c8_file_fmt
      write(6,*)' Filename Format = ',fname_fmt

      do i = 1,i_nbr_files_ret
          call make_fnam_lp(i4times(i),a9_time,istatus)
          filename_in = c_fnames(i)
!         Test whether we want the NetCDF file for this time
          i4time_file_earliest = i4time_sys - i4_acars_window
     1                                      - lag_time_report
          i4time_file_latest =   i4time_sys + i4_acars_window
          
          if(i4times(i) .lt. i4time_file_earliest)then
              write(6,*)
              write(6,*)' File is too early ',a9_time,i
          elseif(i4times(i) .gt. i4time_file_latest)then
              write(6,*)
              write(6,*)' File is too late ',a9_time,i
          else
              write(6,*)
              write(6,*)' File is in time window ',a9_time,i
              write(6,*)' Input file ',filename_in

              if(c8_file_fmt .eq. 'NIMBUS')then ! NIMBUS NetCDF
!                 Read from the ACARS file 
!                 Write to the opened PIN file
!                 ingest_acars_sub.f
                  call get_acars_data(i4time_sys,i4_acars_window
     1                                      ,NX_L,NY_L
     1                                      ,c8_file_fmt
     1                                      ,ext
     1                                      ,l_use_tamdar
     1                                      ,filename_in,istatus)
              elseif(c8_file_fmt .eq. 'WFO')then ! WFO NetCDF
!                 Read from the ACARS file 
!                 Write to the opened PIN file
                  filename13 = cvt_i4time_wfo_fname13(i4times(i))
                  call get_acars_data(i4time_sys,i4_acars_window
     1                                      ,NX_L,NY_L
     1                                      ,c8_file_fmt
     1                                      ,ext
     1                                      ,l_use_tamdar
     1                                      ,filename_in,status)
              elseif(c8_file_fmt .eq. 'AFWA')then ! AFWA ASCII
!                 Read from the ACARS file 
!                 Write to the opened PIN file
                  call get_acars_afwa(i4time_sys,i4_acars_window
     1                                      ,NX_L,NY_L
     1                                      ,ext
     1                                      ,filename_in,istatus)
              else
                  write(6,*)' ERROR, invalid c8_file_fmt: ',c8_file_fmt       
                  istatus = 0
                  return
              endif

          endif
      enddo

!990  close(11) ! Output PIN file

      go to 999

 997  write(6,*)' Error in ACARS ingest (ob time window)'

 999  continue

      write(6,*)' End of acars ingest routine'
 
      return
      end
 
