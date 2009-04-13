
      subroutine ingest_raob(path_to_raw_raob,c8_raob_format,i4time_sys
     1                      ,l_fill_ht,lun_out)

!     Steve Albers FSL   1999       Original Version

!     Input file 
      character*200 filename_in
      character*9 a9_time
      character*180 dir_in
      character*255 c_filespec
      integer max_files
      parameter(max_files = 3000)
      character*255 c_fnames(max_files)
      integer i4times(max_files)
      character*9 a8_to_a9
      character*8 a8_time,a8_time_orig(max_files)
      character*8 c8_raob_format

!     Output file
      character*13 filename13, cvt_i4time_wfo_fname13
      character*31    ext
      integer       len_dir

      character*40 c_vars_req
      character*(*) path_to_raw_raob

      logical l_parse, l_fill_ht

!     Define interval to be used (between timestamps) for creation of SND files
      integer i4_snd_interval
      parameter (i4_snd_interval = 3600)

      iopen = 0

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
          go to 999
      endif

      dir_in = path_to_raw_raob

      call s_len(dir_in,len_dir_in)

      if(c8_raob_format(1:6) .eq. 'NIMBUS')then
          c_filespec = dir_in(1:len_dir_in)//'/*0300o'
          call get_file_times(c_filespec,max_files,c_fnames
     1                       ,i4times,i_nbr_files_ret,istatus)

      elseif(c8_raob_format(1:3) .eq. 'WFO' .or.
     1       c8_raob_format(1:3) .eq. 'RSA'      )then
          c_filespec = dir_in(1:len_dir_in)
          call get_file_times(c_filespec,max_files,c_fnames
     1                       ,i4times,i_nbr_files_ret,istatus)

      elseif(      l_parse(c8_raob_format,'AFGWC')
     1        .or. l_parse(c8_raob_format,'AFWA')   )then
          c_filespec = dir_in(1:len_dir_in)//'/raob.*'
          call get_file_names(c_filespec,i_nbr_files_ret,c_fnames
     1                      ,max_files,istatus)

!         Obtain file times from file names
          do i = 1,i_nbr_files_ret
              call s_len(c_fnames(i),len_fname)
              call get_directory_length(c_fnames(i),len_dir)
              a8_time = c_fnames(i)(len_dir+6:len_fname)
              a8_time_orig(i) = a8_time

              a9_time = a8_to_a9(a8_time)
              call i4time_fname_lp(a9_time,i4times(i),istatus)
              write(6,*)c_fnames(i)(1:len_fname),i4times(i)
          enddo ! i

      elseif(c8_raob_format(1:3) .eq. 'CWB')then
          c_filespec = dir_in(1:len_dir_in)//'/temp*'
          call get_file_names(c_filespec,i_nbr_files_ret,c_fnames
     1                      ,max_files,istatus)

!         Obtain file times from file names
          do i = 1,i_nbr_files_ret
              call s_len(c_fnames(i),len_fname)
              call get_directory_length(c_fnames(i),len_dir)
              a8_time = c_fnames(i)(len_dir+5:len_fname)
              a8_time_orig(i) = a8_time

              a9_time = a8_to_a9(a8_time)
              call i4time_fname_lp(a9_time,i4times(i),istatus)
              write(6,*)c_fnames(i)(1:len_fname),i4times(i)
          enddo ! i

      else
          write(6,*)' Error - Invalid c8_raob_format ',c8_raob_format       
          istatus = 0
          goto 999

      endif

!     Get RAOB Time Window
      call get_windob_time_window('RAOB',i4_wind_ob,istatus)
      if(istatus .ne. 1)goto 997

      call get_tempob_time_window('RAOB',i4_temp_ob,istatus)
      if(istatus .ne. 1)goto 997

      i4_raob_window = max(i4_wind_ob,i4_temp_ob)

!     Loop through raob files and choose ones in time window
      write(6,*)' # of files using filename format ',c8_raob_format      
     1                                        ,' = ',i_nbr_files_ret
      do i = 1,i_nbr_files_ret
          call make_fnam_lp(i4times(i),a9_time,istatus)

          if(c8_raob_format(1:6) .eq. 'NIMBUS')then
              filename_in = dir_in(1:len_dir_in)//'/'//a9_time//'0300o'       
!             i4_raob_window = 60000  ! Temporary for testing
              i4_contains_early = 10800
              i4_contains_late  = 3600

          elseif(c8_raob_format(1:3) .eq. 'WFO')then
              filename13 = cvt_i4time_wfo_fname13(i4times(i))
              filename_in = dir_in(1:len_dir_in)//'/'//filename13      
              i4_contains_early = 43200
              i4_contains_late  = 0

          elseif(c8_raob_format(1:3) .eq. 'RSA')then
              filename13 = cvt_i4time_wfo_fname13(i4times(i))
              filename_in = dir_in(1:len_dir_in)//'/'//filename13      
              i4_contains_early = 43200
              i4_contains_late  = 0

          elseif(      l_parse(c8_raob_format,'AFGWC')
     1            .or. l_parse(c8_raob_format,'AFWA')   )then
              filename_in = dir_in(1:len_dir_in)//'/raob.'//
     1                                            a8_time_orig(i)
              i4_contains_early = 3600
              i4_contains_late  = 0

          elseif(c8_raob_format(1:3) .eq. 'CWB')then 
              filename_in = dir_in(1:len_dir_in)//'/temp'//
     1                      a8_time_orig(i)//'.dat'
              i4_contains_early = 19800         
              i4_contains_late  = 23400       

          else
              write(6,*)' Error - Invalid c8_raob_format '
     1                 ,c8_raob_format    
              istatus = 0
              goto 999

          endif

!         filename_in = 'test.nc                                 '

!         Define limits of RAOB release times we are interested in
          i4time_raob_latest =   i4time_sys + i4_raob_window 
          i4time_raob_earliest = i4time_sys - i4_raob_window 

!         Define limits of file times we are interested in. 
          i4_filetime_latest =   i4time_raob_latest + i4_contains_early
          i4_filetime_earliest = i4time_raob_earliest - i4_contains_late
          
          if(i .eq. 1)then
              write(6,*)' i4 raob sys/window'
     1                   ,i4time_sys,i4_raob_window
              write(6,*)' i4 raob range     '
     1                   ,i4time_raob_earliest,i4time_raob_latest
              write(6,*)' i4 file range     '
     1                   ,i4_filetime_earliest,i4_filetime_latest
          endif

          if(i4times(i) .lt. i4_filetime_earliest)then
              write(6,*)' File is too early ',a9_time,i

          elseif(i4times(i) .gt. i4_filetime_latest)then
              write(6,*)' File is too late ',a9_time,i

          else
              write(6,*)
              write(6,*)' File is in time window ',a9_time,i
              write(6,*)' Input file ',filename_in
              write(6,*)

!             Read from the raw file and write to the opened SND file
              if(c8_raob_format(1:6) .eq. 'NIMBUS' .or.
     1           c8_raob_format(1:3) .eq. 'WFO'         )then

                  call get_raob_data   (i4time_sys,ilaps_cycle_time
     1                ,NX_L,NY_L
     1                ,i4time_raob_earliest,i4time_raob_latest
     1                ,filename_in,lun_out,l_fill_ht,istatus)

              elseif(c8_raob_format(1:3) .eq. 'RSA')then
                  write(6,*)
     1                ' Calling get_rtamps_data (under construction)'       

                  call get_rtamps_data(i4time_sys,ilaps_cycle_time       
     1                ,NX_L,NY_L
     1                ,i4time_raob_earliest,i4time_raob_latest
     1                ,filename_in,lun_out,istatus)

              elseif(      l_parse(c8_raob_format,'AFGWC')
     1                .or. l_parse(c8_raob_format,'AFWA')   )then

!                 Open output SND file 
                  call open_ext(lun_out,i4time_sys,'snd',istatus)

                  call get_raob_data_af(i4time_sys,ilaps_cycle_time
     1                ,NX_L,NY_L
     1                ,i4time_raob_earliest,i4time_raob_latest,a9_time       
     1                ,filename_in,istatus)

              elseif(c8_raob_format(1:3) .eq. 'CWB')then
!                 Open output SND file 
                  call open_ext(lun_out,i4time_sys,'snd',istatus)

                  call get_raob_data_cwb(i4time_sys,ilaps_cycle_time
     1                ,NX_L,NY_L
     1                ,i4time_raob_earliest,i4time_raob_latest,a9_time       
     1                ,filename_in,istatus)

              else
                  write(6,*)' Error - Invalid c8_raob_format '
     1                     ,c8_raob_format
                  istatus = 0
                  goto 999

              endif

          endif

      enddo

      go to 999

 997  write(6,*)' Error in RAOB ingest (ob time windows)'

      go to 999

 998  write(6,*)' Error opening output sounding file: '

 999  continue

      write(6,*)' End of raob ingest routine'

      return
      end
