cdis   
cdis    Open Source License/Disclaimer, Forecast Systems Laboratory
cdis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
cdis    
cdis    This software is distributed under the Open Source Definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    In particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - Redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - Redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - All modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - If significant modifications or enhancements are made to this
cdis    software, the FSL Software Policy Manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
cdis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
cdis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
cdis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
cdis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
cdis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
cdis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
cdis   
cdis
cdis
cdis   
cdis

      subroutine ingest_vad(istatus)

!     Steve Albers      Dec-1997        

!     Input file 
      character*170 filename_in
      character*9 a9_time
      character*150 dir_in
      character*255 c_filespec
      integer max_files
      parameter(max_files = 1000)
      character*255 c_fnames(max_files)
      integer i4times(max_files)

!     Output file
      character*31    ext

      call get_systime(i4time_sys,a9_time,istatus)
      if(istatus .ne. 1)go to 999

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
          istatus = 0
          return
      endif

      iopen = 0

      lag_time_report = 4 * 3600

!     Get List of input /public NetCDF files
      call get_vad_parms(dir_in,istatus)
      if(istatus .eq. 1)then
          write(6,*)'dir_in = ',dir_in
      else
          write(6,*)' Error getting dir_in'
          istatus = 0
          return
      endif

      call s_len(dir_in,len_dir_in)
      c_filespec = dir_in(1:len_dir_in)//'/'//'*0030r'
      call get_file_times(c_filespec,max_files,c_fnames
     1                      ,i4times,i_nbr_files_ret,istatus)

!     Loop through /public NetCDF files and choose ones in time window
      write(6,*)' # of files on /public = ',i_nbr_files_ret
      do i = 1,i_nbr_files_ret
          call make_fnam_lp(i4times(i),a9_time,istatus)
          filename_in = dir_in(1:len_dir_in)//'/'//a9_time//'0030r'

!         Test whether we want the NetCDF file for this time
          i4time_file_earliest = i4time_sys - (ilaps_cycle_time / 2)
          i4time_file_latest =   i4time_sys + (ilaps_cycle_time / 2) 
     1                                      + lag_time_report
          
          if(i4times(i) .lt. i4time_file_earliest)then
              write(6,*)' File is too early ',a9_time,' ',i
          elseif(i4times(i) .gt. i4time_file_latest)then
              write(6,*)' File is too late ',a9_time,' ',i
          else
              write(6,*)' File is in time window ',a9_time,i
              write(6,*)' Input file ',filename_in

              if(iopen .eq. 0)then
!                 Open output PRO file
                  iopen = 1
                  ext = 'pro' ! 'vad'
                  call open_ext(1,i4time_sys,ext(1:3)
     1                                         ,istatus)
                  if(istatus .ne. 1)return
              endif

!             Read from the NetCDF pirep file and write to the opened PIN file
              call get_vad_data(i4time_sys,ilaps_cycle_time
     1                                      ,NX_L,NY_L
     1                                      ,filename_in,istatus)
          endif
      enddo

 999  continue

      write(6,*)' End of ingest_vad subroutine'

      istatus = 1
      return
      end
 
 
       subroutine get_vad_parms(path_to_vad,istatus)

       character*150 path_to_vad
       namelist /vad_nl/ path_to_vad
 
       character*150 static_dir,filename
 
       call get_directory('static',static_dir,len_dir)

       filename = static_dir(1:len_dir)//'/vad.nl'
 
       open(21,file=filename,status='old',err=900)
       read(21,vad_nl,err=901)
       close(21)

       istatus = 1
       return

  900  print*,'error opening file ',filename
       istatus = 0
       return

  901  print*,'error reading vad_nl in ',filename
       istatus = 0
       return

       end
