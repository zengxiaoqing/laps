
      subroutine ingest_pireps(istatus)

!     Steve Albers                      Original Version
!     Ken Dritz     15-Jul-1997         Added call to get_grid_dim_xy to
!                                       get the values of NX_L, NY_L
!     Ken Dritz     15-Jul-1997         Pass NX_L, NY_L to get_pirep_data

!     Input file 
      character*170 filename_in
      character*9 a9_time
      character*180 dir_in
      character*255 c_filespec
      integer max_files
      parameter(max_files = 3000)
      character*255 c_fnames(max_files)
      integer i4times(max_files)
      character*8 c8_project

!     Output file
      character*13 filename13
      character*31    ext
      integer       len_dir

      character*40 c_vars_req
      character*180 c_values_req

      call get_systime(i4time,a9_time,istatus)
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

      

      call get_laps_cycle_time(ilaps_cycle_time,istatus)
      if(istatus .eq. 1)then
          write(6,*)' ilaps_cycle_time = ',ilaps_cycle_time
      else
          write(6,*)' Error getting laps_cycle_time'
          return
      endif

      lag_time_report = 3600

!     Get List of input /public NetCDF files
      c_vars_req = 'path_to_raw_pirep'
      call get_static_info(c_vars_req,c_values_req,1,istatus)
      if(istatus .eq. 1)then
          write(6,*)c_vars_req(1:30),' = ',c_values_req
          dir_in = c_values_req
      else
          write(6,*)' Error getting ',c_vars_req
          return
      endif

!     dir_in = path_to_raw_pirep

      call s_len(dir_in,len_dir_in)


!     if(c8_project .eq. 'WFO')then
        c_filespec = dir_in(1:len_dir_in)
!     else
!       c_filespec = dir_in(1:len_dir_in)//'*0005r'
!     endif

      call get_file_times(c_filespec,max_files,c_fnames
     1                      ,i4times,i_nbr_files_ret,istatus)

!     Open output PIN file
      if(i_nbr_files_ret .gt. 0 .and. istatus .eq. 1)then
          ext = 'pin'
          if(istatus .ne. 1)then
              write(6,*)' Error opening output file'
              go to 999
          endif
      
      else
          write(6,*)' No raw data files identified'
          go to 999
      endif


!     Loop through /public NetCDF files and choose ones in time window
      write(6,*)' # of files = ',i_nbr_files_ret
      do i = 1,i_nbr_files_ret
          call make_fnam_lp(i4times(i),a9_time,istatus)
          if(c8_project .eq. 'WFO')then
            filename_in = c_fnames(i)
          else
            filename_in = dir_in(1:len_dir_in)//a9_time//'0005r'
          endif
!         filename_in = 'test.nc                                 '

!         Test whether we want the NetCDF file for this time
          i4time_file_latest = i4time + (ilaps_cycle_time / 2) 
     1                                + lag_time_report
          i4time_file_earliest = i4time - (ilaps_cycle_time / 2)
          
          if(i4times(i) .lt. i4time_file_earliest)then
              write(6,*)' File is too early ',a9_time,i
          elseif(i4times(i) .gt. i4time_file_latest)then
              write(6,*)' File is too late ',a9_time,i
          else
              write(6,*)' File is in time window ',a9_time,i
              write(6,*)' Input file ',filename_in

!             Read from the NetCDF pirep file and write to the opened PIN file
              if(c8_project .eq. 'WFO')then
                write(6,*)' Beginning AWIPS PIREP ingest - ',filename_in
                write(6,*)
                call get_pirep_data_WFO(i4time,ilaps_cycle_time
     1                                      ,filename_in
     1                                      ,ext
     1                                      ,NX_L,NY_L
     1                                      ,istatus)
                write(6,*)
                write(6,*)' End of AWIPS PIREP ingest - ', filename_in
                write(6,*)
              else
                call get_pirep_data(i4time,ilaps_cycle_time,filename_in
     1                                      ,ext
     1                                      ,NX_L,NY_L
     1                                      ,istatus)
              endif
          endif
      enddo


!     close(11) ! Output PIN file

 999  continue

      write(6,*)' End of pirep ingest routine'

      return
      end
 
