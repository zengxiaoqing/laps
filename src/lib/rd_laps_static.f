      subroutine rd_laps_static (dir,laps_dom_file,imax,jmax,n_grids,
     1                           var,units,comment,data,grid_spacing,
     1                           status)

      implicit none

C This routine calls functions in static_routines.c.  To create the new
C static file, a cdl file named <laps_dom_file>.cdl (ie. nest7grid.cdl)
C must be located in the LAPS cdl directory and is written out to the
C directory specified in "dir".

C Passed in variables
      character      dir*(*), laps_dom_file*(*)
      integer        imax, jmax, n_grids
      character*(*)  var(n_grids)
      character*(*)  units(n_grids)
      character*(*)  comment(n_grids)
      real         data(imax,jmax,n_grids), grid_spacing
      integer        status

C Local variables

      character*150   file_name
      integer        var_len, com_len, unit_len,
     1               ERROR(3),
     1               no_laps_diag,
     1               flag

      integer        f_len

      logical        l_some_missing

      COMMON         /PRT/FLAG

      ERROR(1)=1
      ERROR(2)=0
      ERROR(3)=-2
      l_some_missing = .false.

C  BEGIN SUBROUTINE

      call make_static_fname(dir,laps_dom_file,file_name,f_len,status)
      if (status .eq. ERROR(2)) goto 980
      
      var_len = len(var(1))
      com_len = len(comment(1))
      unit_len = len(units(1))

      call read_cdf_static(file_name,f_len,var,var_len,comment,com_len,
     1                     units,unit_len,imax,jmax,n_grids,data,
     1                     grid_spacing,no_laps_diag,status)

      if (status .ge. 0) then       !return from read with no errors

c       do i = 1, n_grids
c         comment(i) = c_comment(i)
c         units(i) = c_units(i)
c       enddo

        if (status .gt. 0) l_some_missing = .true.

        if(l_some_missing) then
          STATUS=ERROR(3) != -2
        else
          STATUS=ERROR(1) !=  1
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
     1write(6,*) 'Error in imax,jmax, or n_grids...read aborted.'
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
980   if (flag .ne. 1) then
        write(6,*) 'Length of dir+file-name is greater than 150 char.'
        write(6,*) 'Static file cannot be accessed.'
      endif
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

      integer end_dir, end_dom,
     1          ERROR(2),
     1          status

      integer f_len

      character       dir*(*)       !Directory to read data from
      character       laps_dom_file*(*)
      character*30    dn_laps_dom   !downcase of laps_dom_file

      character*(*)    file_name

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
C ******  find end of file_name
C
      f_len = len(file_name)
C
C ****  make file_name
C
      if (end_dir+end_dom+7 .gt. f_len) then
        status = ERROR(2)
        goto 999
      else
        file_name = dir(1:end_dir)//'static.'//dn_laps_dom(1:end_dom)
        f_len = end_dir+7+end_dom
      endif
C
C ****  Return normally.
C
      status = ERROR(1)
999   return
      end

C##########################################################################
