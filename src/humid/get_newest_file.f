      subroutine get_newest_file (filename,time_diff,
     1     path,pi,file,ext,ei,istatus)

      implicit none

c     parameters
      
      integer maxfiles
      parameter (maxfiles = 3000)
      
c     parameter list variables

      character*9 filename 
      integer time_diff
      character*512 path
      integer pi
      character*9 file
      character*120 ext
      integer ei
      integer istatus

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

      
c     get the  filenames available from get_filenames

      call get_file_names (path(1:pi), numoffiles, filenames,
     1     maxfiles, istatus)


      if (istatus.ne.1) then    !failure, return with failure code
         return
      endif

      maxhash = ' '
      hash = ' '

      do i = 1,numoffiles
         hashpoint = index (filenames(i),'tpw')
         if (hashpoint .ne. 0) hash = filenames(i)
         if (hash.gt.maxhash) maxhash = hash
      enddo                     ! i
   
      hashpoint = index (hash,'.')

      if (hashpoint .eq. 0) then
         istatus = 0
         return
      endif

      hash = maxhash

      file_new = hash(hashpoint-9:hashpoint)
      ext = hash(hashpoint+1:256)

      call s_len(ext, ei, istatus)

c     determine i4 times of pertenant files

      call i4time_fname_lp (file_new, i4time_found, istatus)
      call i4time_fname_lp (filename, i4time_input, istatus)

c     compute time difference
      
      i4time_difference = i4time_input-i4time_found

      if (i4time_difference .lt. time_diff) then !success in finding file
         file = file_new
         write(6,*) 'Found file ', file_new//'.'//ext(1:ei)
         write(6,*) 'Age of GVAP file is ',i4time_difference, ' seconds'
         return
      else
         istatus = 0 
         write(6,*)'GVAP file not found'
         return
      endif

      return
      end
