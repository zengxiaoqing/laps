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
      
      subroutine get_newest_file (filename,time_diff,file_name_length,
     1     path,pi,file,de,dei,ext,ei,istatus)

      implicit none

c     parameters
      
      integer maxfiles
      parameter (maxfiles = 5000)
      
c     parameter list variables

      character*9 filename 
      integer time_diff
      character*256 path
      integer pi
      character*9 file
      character*120 ext
      integer ei
      integer istatus
      integer file_name_length
      character*120 de
      integer dei

c     internal variables
      character*256 filenames (maxfiles)
      character*256 hash, maxhash
      integer hashpoint
      integer numoffiles
      integer i4time_input
      integer i4time_found
      integer i4time_difference
      character*9 file_new
      integer i

c     code

      call s_len(de, dei)

      
c     get the  filenames available from get_filenames
      write (6,*) 'Calling get_file_names with path:'
      write (6,*) path(1:pi)
      call get_file_names (path(1:pi), numoffiles, filenames,
     1     maxfiles, istatus)
      write (6,*) 'Get_file_names returns with:'
      write (6,*) 'numoffiles: ',numoffiles
      write (6,*) 'getfilenames istatus = ',istatus


      if (istatus.ne.1) then    !failure, return with failure code
         return
      endif

      maxhash = ' '
      hash = ' '
      
      do i = numoffiles,1,-1    !,numoffiles
         hashpoint = index (filenames(i),de(1:dei))
         if (hashpoint .ne. 0) hash = filenames(i)
         hashpoint = index (hash,'.')
         if (hashpoint .eq. 0) go to 15 ! try next file
         file_new = hash(hashpoint-file_name_length:hashpoint)
         ext = hash(hashpoint+1:256)
         call s_len(ext, ei)
      
c     determine i4 times of pertenant files
      
         call i4time_fname_lp (file_new, i4time_found, istatus)
         call i4time_fname_lp (filename, i4time_input, istatus)
      
c     compute time difference
      
         i4time_difference = i4time_input-i4time_found
         i4time_difference  = (i4time_difference)
      
         if (i4time_difference .lt. time_diff 
     1                 .and. i4time_difference .gt. 0) then !success in finding file
            file = file_new
            write(6,*) 'Found file ', file_new//'.'//ext(1:ei)
            write(6,*) 
     1    'Age of file found is ',i4time_difference,' seconds'
            return
         endif
 15      continue
      enddo                     ! i done with all files
      istatus  = 0
      write (6,*) 'File not found in Get_newest_file'
      return
      end
