cdis    Forecast Systems Laboratory
cdis    NOAA/OAR/ERL/FSL
cdis    325 Broadway
cdis    Boulder, CO     80303
cdis
cdis    Forecast Research Division
cdis    Local Analysis and Prediction Branch
cdis    LAPS
cdis
cdis    This software and its documentation are in the public domain and
cdis    are furnished "as is."  The United States government, its
cdis    instrumentalities, officers, employees, and agents make no
cdis    warranty, express or implied, as to the usefulness of the software
cdis    and documentation for any purpose.  They assume no responsibility
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis
cdis    Permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  All modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making
cdis    the modifications.  If significant modifications or enhancements
cdis    are made to this software, the FSL Software Policy Manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis
cdis
cdis
cdis
cdis
cdis
cdis
      subroutine get_file_names(pathname_in,numoffiles,c_filenames,
     1     max_files,istatus)

c
c     Author : Dan Birkenheuer
c
c     Date: 5/22/95
c     12/2/96 modified error output. DB
c     
c     to be used with C-routine getfilenames_c.. this is the Fortran
c     wrapper for that routine.
c

      implicit none

      integer numoffiles,max_files
      character*(*) pathname_in
      character*(*) c_filenames(max_files)
      
      
      character*256 file_names(20000)
      character*256 dirpath
      character*256 filter
      integer istatus,i,nc,nindx,status,nindx2
      
      nindx = index (pathname_in, '*')
      
      if( nindx .eq. 0) then
         
         dirpath = pathname_in
         filter = '*'
         
      else
         
         nc = len(pathname_in)
         i = nc
         do while (pathname_in(i:i) .ne. '/' .and. i .gt. 0)
            i = i-1
            if (i.eq.0) goto 122
         enddo
 122     continue
         
         if (i.gt.0) then
            dirpath = pathname_in(1:i)
            filter = pathname_in(i+1:nc)
         else
            dirpath = '.'
            filter = pathname_in(i+1:nc)
         endif
         
        endif
        
        
        call getfilenames_c (dirpath,file_names,numoffiles,filter,
     1       status)
        
        
        if(numoffiles.gt.max_files) then
           print *, 'Error calling Get_file_names'
           print *, 'number of files returned exceeds array allocation'
           print *, ' for file_names.  (numoffiles > max_files)'
           istatus = 0          ! failure
           print *, 'Maxfile is ',max_files
           print *, 'Number of files returned is ',numoffiles
           return
        endif
        call s_len(dirpath,nindx2)
c     print*,'>',dirpath,'<',
c     print*,'nindx2 ',nindx2,string_space(dirpath)
        
        do i = 1, numoffiles
c     nindx = index(file_names(i), ' ')
           call s_len(file_names(i),nindx)
           if( dirpath(nindx2:nindx2) .ne. '/')then
              c_filenames(i) = dirpath(1:nindx2)//'/'
     1             //file_names(i)(1:nindx)
           else
              c_filenames(i) = dirpath(1:nindx2)
     1             //file_names(i)(1:nindx)
           endif
           
        enddo
        
        istatus = status
        
        if (istatus.ne.1) then
           return
        endif
        
c     sort data
        
        call sort_fn (c_filenames,numoffiles,istatus)
        
        return
        end
      
      
      
      
      
      subroutine sort_fn (names,number,istatus)
      
      integer number            ! Input
      integer istatus           ! Output
      character*(*) names(number) ! Input/Output
      
      integer i4time(number)    ! Local
      character*500 names_buf   ! Local
      logical l_switch          ! Local
      
!       This does a Bubble Sort on the list of names.
!       Note that this is relatively inefficient, more efficient
!       alternatives include a "tree" sort or various C sorting utilities.
!       WARNING: Input names must be <= 500 characters in length.

!       Steve Albers    June 1995       Original Version
!       Steve Albers    July 1998       Fix for Y2K problem

      if(number .eq. 0)then
         istatus = 1
         return
      endif
      
      n_non_numeric = 0

!       Set up array of i4times for the files, 0 is used for non-numeric files

      do i = 1,number

         call i4time_fname_lp(names(i),i4time(i),istatus)
         
         if(istatus .ne. 1)then
            i4time(i) = 0       ! File is non-numeric and has no i4time
            n_non_numeric = n_non_numeric + 1
         endif
         
      enddo                     ! i

 10   l_switch = .false.

      do i=1,number-1
         
!           Compare two names
         if( (i4time(i) .gt. i4time(i+1) ) .OR.
     1        (i4time(i) .eq. i4time(i+1) .and. 
     1        names(i)  .gt. names(i+1)  )      )then
            
!               Switch the names
            names_buf = names(i)
            names(i) = names(i+1)
            names(i+1) = names_buf
            
            i4time_buf = i4time(i)
            i4time(i) = i4time(i+1)
            i4time(i+1) = i4time_buf
            
            l_switch = .true.

         endif

      enddo
      
      if(l_switch)then
!           write(6,*)' Sort_fn: switched one or more names'
         go to 10
      else
         write(6,*)' Sort_fn: names sorted,'
     1        ,' non_numeric/total = ',n_non_numeric,number
      endif
      
      istatus = 1
      
      return
      end

