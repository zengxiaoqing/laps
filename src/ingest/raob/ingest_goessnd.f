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

      subroutine ingest_goessnd (path_to_data, i4time_sys, 
     1     lun_out, istatus)


      implicit none
      
      character*(*) path_to_data
      integer i4time_sys, lun_out, istatus


      integer max_files
      parameter (max_files = 3000)
      integer i4time_begin,i4time_end ! bracket of i4time of interest
      character*256 file_names(max_files) !file names found in path of interest
      character*256 files_to_process(max_files) ! files meeting time criteria
      character*9 a9_test !a9 filename to test for criteria
      integer     i4time_test   ! i4time of filename to test for time criteria
      integer num_files
      integer fname_len         ! length of file name found in path
      integer path_len          ! length of path_to_data
      integer i,j,ii,jj,lct
      real mdf               !missing data flag
      character*50 ext
      character *256 directory
      integer len_dir

      istatus = 0

      

c     get missing data flag
      call get_r_missing_data (mdf,istatus) 

c     define the spatial window for data
c     obtain dimensions of laps domain ii,jj
      call get_grid_dim_xy(ii,jj,istatus)


c     define time window to look for data
      i4time_begin = i4time_sys -4200
      i4time_end   = i4time_sys +3600

c     get the filenames to consider processing
      call s_len (path_to_data,path_len)
      write (6,*) path_to_data(1:path_len)
      call get_file_names (path_to_data(1:path_len), num_files,
     1     file_names,max_files,istatus)

      if(istatus.eq.1) then
         
c     trim file names removing path
         do i=1,num_files
            call s_len(file_names(i),fname_len)
            file_names(i) = file_names(i)(path_len+1:fname_len)
         enddo

      else
         write(6,*) 'error in get_file_names'
         return
      endif

      
c     go through list of filenames and test their valid times against window

      j = 0 ! counter to increment for good files to run through scanner
      do i = 1, num_files
         a9_test = file_names(i)(5:9)//file_names(i)(11:14)
         call i4time_fname_lp(a9_test,i4time_test,istatus)


c         i4time_begin = i4time_begin - 3800 ! inserted to expand time window

         if (i4time_test.ge.i4time_begin .and.
     1        i4time_test.le.i4time_end)then
            j = j+1
            files_to_process(j) = file_names(i)
         endif
      enddo                     !i                    

      if (j.eq.0) then
         write (6,*) 'no files to process'
         istatus = 0
         return
      endif

c     New call for Steve Albers check of lun
      call open_ext ( lun_out, i4time_sys, 'snd', istatus)

c     what is the new file name length (same for all files)

      call s_len(files_to_process(1),fname_len)

      do i = 1,j                !for all files to be processed

         write (6,*) 'processing data GOES data files ',j
         write (6,*) 'filename is ',files_to_process(i)
         call process_goes_snd (path_to_data, path_len, 
     1        files_to_process(i),
     1        fname_len, ii,jj,
     1        i4time_begin,i4time_end,mdf,lun_out,istatus)

      enddo                     ! i

      istatus = 1
      return
      end




