      subroutine rd_laps_static (dir,laps_dom_file,imax,jmax,n_grids,
     1                           var,units,comment,data,grid_spacing,
     1                           status)

      implicit none

C MAX_GRIDS is the max number of grids that may be passed in for writing
C out to the static file.  An additional grid, ZIN, is generated within
C this subroutines from the AVG grid
      integer*4      max_grids
      parameter      (MAX_GRIDS=7)

      character      dir*(*),
     1               laps_dom_file*(*)

      integer*4      imax,
     1               jmax,
     1               n_grids

      character*3    var(n_grids)
      character*4    c_var(8)
      character*125  comment(8)
      character*126  c_comment(8)
      character*10   units(n_grids)
      character*11   c_units(8)
      character*91   file_name

      real*4         data(imax,jmax,n_grids),
     1               grid_spacing

      integer*4      status

C Local variables

      integer*4      i,
     1               ERROR(3),
     1               no_laps_diag,
     1               flag

      integer*2      f_len

      logical        l_some_missing

      COMMON         /PRT/FLAG

      ERROR(1)=1
      ERROR(2)=0
      ERROR(3)=-2
      l_some_missing = .false.

C  BEGIN SUBROUTINE

      do i=1,n_grids
        call upcase(var(i),c_var(i))
      enddo

      call make_static_fname(dir,laps_dom_file,file_name,f_len,status)

      print *,'rd_laps_static: ',file_name

      call read_cdf_static(file_name,f_len,c_var,c_comment,c_units,
     1                     imax,jmax,n_grids,data,grid_spacing,
     1                     no_laps_diag,status)

      if (status .ge. 0) then       !return from read with no errors

        do i = 1, n_grids
          comment(i) = c_comment(i)
          units(i) = c_units(i)
        enddo

        if (status .gt. 0) l_some_missing = .true.

        if(l_some_missing) then
          STATUS=ERROR(3)
        else
          STATUS=ERROR(1)
        endif
      endif

      if (status .eq. -1) GOTO 940 !error opening file
      if (status .eq. -3) GOTO 950 !error in imax,jmax or n_grids
      if (status .eq. -4) GOTO 960 !error reading file data
      if (status .eq. -5) GOTO 970 !error reading header
C
C ****  Return normally.
C
      status = ERROR(1)
999   return
C
C ****  Error trapping.
C
940   if (flag .ne. 1)
     1   write(6,*) 'Error opening file to be read.'
      status = ERROR(2)
      goto 999
C
950   if (flag .ne. 1)
     1write(6,*) 'Error in imax,jmax, or n_grids...write aborted.'
      status = ERROR(2)
      goto 999
C
960   if (flag .ne. 1)
     1   write(6,*) 'Error reading file data.'
      status = ERROR(2)
      goto 999
C
970   if (flag .ne. 1)
     1   write(6,*) 'Error reading header info.'
      status = ERROR(2)
      goto 999
C
      end

C########################################################################
      subroutine make_static_fname(dir,laps_dom_file,file_name,
     1                             f_len,status)
C
C**********************************************************************
C
C      Subroutine MAKE_STATIC_FNAME
C
C      Author:    Linda Wharton 4/93
C
C      Inputed DIR and LAPS_DOM_FILE are converted to ASCII values and
C      B_FILENAME is created.
C
C**********************************************************************

      IMPLICIT  NONE

      integer*4 end_dir, end_dom,
     1          ERROR(2),
     1          i, j,
     1          status

      integer*2 f_len,length

      character       dir*(*)       !Directory to read data from
      character       laps_dom_file*(*)
      character*30    dn_laps_dom   !downcase of laps_dom_file

      Logical         space

      character*91    file_name

      ERROR(1)=1
      ERROR(2)=0

      call downcase(laps_dom_file,dn_laps_dom)

C ******  find end of dir
C
      call s_len(dir,end_dir)
C
C ******  find end of laps_dom_file
C
      call s_len(laps_dom_file,end_dom)
C
C ****  make file_name
C
      file_name = dir(1:end_dir)//'static.'//dn_laps_dom(1:end_dom)
      f_len = end_dir+7+end_dom
C
C ****  Return normally.
C
      status = ERROR(1)
999   return
      end

C##########################################################################
