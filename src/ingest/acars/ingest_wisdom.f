
      subroutine ingest_wisdom(istatus)

      use mem_namelist, ONLY: path_to_wisdom

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

!     Output file
      character*31    ext

      character*40 c_vars_req
      character*180 c_values_req

      lun_out = 31

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

      lag_time_report = 3600

      call s_len(path_to_wisdom,lenp)
      write(6,*)' path_to_wisdom = ',path_to_wisdom(1:lenp)
      dir_in = path_to_wisdom

      ext_in = 'wfo'
      c8_file_fmt = 'WFO'
      call s_len(dir_in,len_dir_in)
      c_filespec = dir_in(1:len_dir_in)

      fname_fmt = ext_in

!     Get list of files
      call get_file_times(c_filespec,max_files,c_fnames
     1                      ,i4times,i_nbr_files_ret,istatus)

!     Open output PIN file for appending
      if(i_nbr_files_ret .gt. 0 .and. istatus .eq. 1)then
          ext = 'pin'
      else
          write(6,*)' No raw data files identified of fname fmt '
     1              ,fname_fmt
          goto999
      endif

!     Get WISDOM Time Window
      call get_windob_time_window('ACARS',i4_wind_ob,istatus)
      if(istatus .ne. 1)goto 997

      i4_wisdom_window = i4_wind_ob
      i4time_earliest = i4time_sys - i4_wisdom_window
      i4time_latest   = i4time_sys + i4_wisdom_window

!     Loop through WISDOM files and choose ones in time window
      write(6,*)' # of raw WISDOM files = ',i_nbr_files_ret
      write(6,*)' File Format = ',c8_file_fmt
      write(6,*)' Filename Format = ',fname_fmt

      do i = 1,i_nbr_files_ret
          call make_fnam_lp(i4times(i),a9_time,istatus)
          filename_in = c_fnames(i)

!         Test whether we want the NetCDF file for this time
          i4time_file_earliest = i4time_sys - i4_wisdom_window
     1                                      - lag_time_report
          i4time_file_latest =   i4time_sys + i4_wisdom_window
          
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

!             Read from the WISDOM file 
!             Write to the opened PIN file
              call get_wisdom_data(i4time_sys,ilaps_cycle_time        ! I
     1                                      ,NX_L,NY_L                ! I
     1                                      ,i4time_earliest          ! I
     1                                      ,i4time_latest            ! I
!    1                                      ,c8_file_fmt
!    1                                      ,ext
     1                                      ,filename_in              ! I
     1                                      ,lun_out                  ! I
     1                                      ,istatus)                 ! I
              if(istatus .ne. 1)then
                  write(6,*)' bad istatus from get_wisdom_data ',istatus      
              endif
          endif
      enddo

      go to 999

 997  write(6,*)' Error in WISDOM ingest (ob time window)'

 999  continue

      write(6,*)' End of wisdom ingest routine'
 
      return
      end
 
