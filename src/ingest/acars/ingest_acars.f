
!     Ken Dritz      15-Jul-1997        Added call to get_grid_dim_xy to get
!                                       the values of NX_L, NY_L
!     Ken Dritz      15-Jul-1997        Pass NX_L, NY_L to get_acars_data
!     Steve Albers      Dec-1997        Call to 'open_lapsprd_file_append'

!     Input file 
      character*70 filename_in
      character*9 a9_time
      character*80 dir_in
      character*255 c_filespec
      integer max_files
      parameter(max_files = 3000)
      character*255 c_fnames(max_files)
      integer i4times(max_files)

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

      lag_time_report = 3600

!     Open output PIN file
      ext = 'pin'
      call open_lapsprd_file_append(11,i4time_sys,ext(1:3),istatus)

!     Get List of input /public NetCDF files
      c_vars_req = 'path_to_qc_acars'
      call get_static_info(c_vars_req,c_values_req,1,istatus)
      if(istatus .eq. 1)then
          write(6,*)c_vars_req(1:30),' = ',c_values_req
          dir_in = c_values_req
      else
          write(6,*)' Error getting ',c_vars_req
          stop   
      endif

      call s_len(dir_in,len_dir_in)
      c_filespec = dir_in(1:len_dir_in)//'*0q.cdf'
      call get_file_times(c_filespec,max_files,c_fnames
     1                      ,i4times,i_nbr_files_ret,istatus)

!     Loop through /public NetCDF files and choose ones in time window
      write(6,*)' # of files on /public = ',i_nbr_files_ret
      do i = 1,i_nbr_files_ret
          call make_fnam_lp(i4times(i),a9_time,istatus)
          filename_in = dir_in(1:len_dir_in)//a9_time//'q.cdf'
!         filename_in = 'test.nc                                 '

!         Test whether we want the NetCDF file for this time
          i4time_file_earliest = i4time_sys - (ilaps_cycle_time / 2)
     1                                      - lag_time_report
          i4time_file_latest =   i4time_sys + (ilaps_cycle_time / 2) 
          
          if(i4times(i) .lt. i4time_file_earliest)then
              write(6,*)' File is too early ',a9_time,i
          elseif(i4times(i) .gt. i4time_file_latest)then
              write(6,*)' File is too late ',a9_time,i
          else
              write(6,*)' File is in time window ',a9_time,i
              write(6,*)' Input file ',filename_in

!             Read from the NetCDF pirep file and write to the opened PIN file
              call get_acars_data(i4time_sys,ilaps_cycle_time
     1                                      ,NX_L,NY_L
     1                                      ,filename_in,istatus)
          endif
      enddo

 990  close(11) ! Output PIN file

 999  continue

      write(6,*)' End of acars ingest program'

      end
 
