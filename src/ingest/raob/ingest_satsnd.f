
      subroutine ingest_satsnd(path_to_raw_satsnd)

!     Steve Albers      May-1999       Original Version

!     Input file 
      character*170 filename_in
      character*9 a9_time
      character*180 dir_in
      character*6 ext_in
      character*255 c_filespec
      integer max_files
      parameter(max_files = 3000)
      character*255 c_fnames(max_files)
      integer i4times(max_files)

!     Output file
      character*31    ext

      character*(*) path_to_raw_satsnd

      call get_systime(i4time_sys,a9_time,istatus)
      if(istatus .ne. 1)go to 999

!     i4time = (i4time/3600) * 3600

      call get_grid_dim_xy(NX_L,NY_L,istatus)
      if (istatus .ne. 1) then
          write (6,*) 'Error getting horizontal domain dimensions'
          go to 999
      endif
 
      lag_time_report = 3600

      dir_in = path_to_raw_satsnd

      ext_in = 'satsnd'

      call s_len(dir_in,len_dir_in)
      c_filespec = dir_in(1:len_dir_in)//'/'//'*00.'//ext_in

!     Get list of files
      call get_file_times(c_filespec,max_files,c_fnames
     1                      ,i4times,i_nbr_files_ret,istatus)

      call get_tempob_time_window('SATSND',i4_satsnd_window,istatus)
      if(istatus .ne. 1)return

      lun_in = 21
      lun_out = 11

!     Open output SND file for appending
      if(i_nbr_files_ret .gt. 0 .and. istatus .eq. 1)then
          ext = 'snd'
          call open_lapsprd_file_append(lun_out,i4time_sys,ext(1:3)
     1                                                    ,istatus)       
      else
          write(6,*)' No raw data files identified:',' *.',ext_in
          goto999
      endif

!     Loop through SATSND files and choose ones in time window
      write(6,*)' # of raw SATSND files = ',i_nbr_files_ret
     1         ,' *.',ext_in
      do i = 1,i_nbr_files_ret
          call make_fnam_lp(i4times(i),a9_time,istatus)
          filename_in = c_fnames(i)

!         Test whether we want the NetCDF file for this time
          i4time_file_earliest = i4time_sys - i4_satsnd_window
     1                                      - lag_time_report
          i4time_file_latest =   i4time_sys + i4_satsnd_window
          
          if(i4times(i) .lt. i4time_file_earliest)then
              write(6,*)' File is too early ',a9_time,i
          elseif(i4times(i) .gt. i4time_file_latest)then
              write(6,*)' File is too late ',a9_time,i
          else
              write(6,*)' File is in time window ',a9_time,i
              write(6,*)' Input file ',filename_in

!             Read from the SATSND file 
!             Write to the opened SND file
              call get_satsnd_afwa(i4time_sys,i4_satsnd_window
     1                            ,NX_L,NY_L
     1                            ,lun_in,filename_in,lun_out,istatus)       
          endif
      enddo

 990  close(lun_out) ! Output SND file

 999  continue

      write(6,*)' End of satsnd ingest routine'

      return
      end
 
