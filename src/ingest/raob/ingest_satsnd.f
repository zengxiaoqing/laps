
      subroutine ingest_satsnd(path_to_raw_satsnd,c8_project
     1                        ,i4time_sys,lun_out)

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
      character*8 c8_project

!     i4time = (i4time/3600) * 3600

      call get_grid_dim_xy(NX_L,NY_L,istatus)
      if (istatus .ne. 1) then
          write (6,*) 'Error getting horizontal domain dimensions'
          go to 999
      endif

      dir_in = path_to_raw_satsnd

      write(6,*)' Subroutine ingest_satsnd - path is: ',trim(dir_in)

      ext_in = 'satsnd'

      call s_len(dir_in,len_dir_in)

      if(c8_project .eq. 'AFWA')then
          if(len_dir_in .gt. 0)then
              c_filespec = dir_in(1:len_dir_in)//'/'//'*00.'//ext_in
          else
              c_filespec = ' '
          endif

          i4_contains_early = 0
          i4_contains_late  = 0
          lag_time_report = 3600

      else ! MADIS POES
          if(len_dir_in .gt. 0)then
              c_filespec = dir_in(1:len_dir_in)//'/'//'*00'
          else
              c_filespec = ' '
          endif

          i4_contains_early = 0
          i4_contains_late  = 3600
          lag_time_report = 0 ! Pretty much an obsolete variable?

      endif

!     Get list of files
      call get_file_times(c_filespec,max_files,c_fnames
     1                   ,i4times,i_nbr_files_ret,istatus)

      call get_tempob_time_window('SATSND',i4_satsnd_window,istatus)
      if(istatus .ne. 1)return

      lun_in = 22

!     Open output SND file for appending
      if(i_nbr_files_ret .gt. 0 .and. istatus .eq. 1)then
          ext = 'snd'
          call open_ext(lun_out,i4time_sys,ext(1:3),istatus)
      else
          write(6,*)' No raw data files identified:',' *.',ext_in
          write(6,*)' c_filespec = ',trim(c_filespec)
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
     1                                      - i4_contains_late     
     1                                      - lag_time_report

          i4time_file_latest =   i4time_sys + i4_satsnd_window
     1                                      + i4_contains_early     
          
          if(i4times(i) .lt. i4time_file_earliest)then
              write(6,*)' File is too early ',a9_time,i
          elseif(i4times(i) .gt. i4time_file_latest)then
              write(6,*)' File is too late ',a9_time,i
          else
              write(6,*)' File is in time window ',a9_time,i
              write(6,*)' Input file ',filename_in

!             Read from the SATSND file 
!             Write to the opened SND file
              if(c8_project .eq. 'AFWA')then
                  call get_satsnd_afwa(i4time_sys,i4_satsnd_window
     1                                ,NX_L,NY_L,lun_in,filename_in
     1                                ,lun_out,istatus)       
              else
                  call get_poes_data
     +                   (i4time_sys,ilaps_cycle_time,NX_L,NY_L
     +                   ,i4time_earliest,i4time_latest
     +                   ,filename_in
     +                   ,lun_out
     +                   ,istatus)
              endif
          endif
      enddo

 990  continue

!     close(lun_out) ! Output SND file

 999  continue

      write(6,*)' End of satsnd ingest routine'

      return
      end
 
