      subroutine wrt_laps_static (dir,laps_dom_file,imax,jmax,n_grids,
     1                            var,comment,data,zin,model,
     1                            grid_spacing,status)

      implicit none

      include 'lapsparms.for'

C MAX_GRIDS is the max number of grids that may be passed in for writing
C out to the static file.  An additional grid, ZIN, is generated within
C this subroutines from the AVG grid

C This routine calls functions in static_routines.c.  To create the new
C static file, a cdl file named <laps_dom_file>.cdl (ie. nest7grid.cdl)
C must be located in the same directory as the output file will be written
C to.  Generally, this is in laps/nest7grid/static.

      integer*4      max_grids
      parameter      (MAX_GRIDS=7)

      character      dir*(*),
     1               laps_dom_file*(DOMAIN_NAME_LEN)

      integer*4      imax,
     1               jmax,
     1               n_grids

      character*3    var(MAX_GRIDS)
      character*125  comment(MAX_GRIDS)
      character*92   file_name

      real*4         data(imax,jmax,n_grids),
     1               zin(imax,jmax),
     1               ZtoPsa,
     1               grid_spacing

      character*131  model

      integer*4      i4time_now_gg,
     1               status

C Local variables

      character*3    var_in(MAX_GRIDS)
      character*24   asctime

      integer*4      i4time,
     1               i, iz, jz, z,
     1               ERROR(2),
     1               flag

      integer*2      f_len

      real*4         psa

      COMMON         /PRT/flag

      ERROR(1)=1
      ERROR(2)=0

C  BEGIN SUBROUTINE

      if (n_grids .gt. MAX_GRIDS) goto 930

      do i=1,n_grids
        call upcase(var(i),var_in(i))
      enddo

      call make_static_fname(dir,laps_dom_file,file_name,f_len,status)

      i4time = i4time_now_gg()
      call cv_i4tim_asc_lp(i4time,asctime,status)

C generate zin grid
C   find index of AVG
      do i = 1, n_grids
        if (var(i) .eq. 'AVG') z = i
      enddo

      do iz = 1, imax
        do jz= 1, jmax
          psa = ZtoPsa(data(iz,jz,z))
          zin(iz,jz) = (20.0 - ((psa - 100.0) * 0.02))
        enddo
      enddo

      call write_cdf_static(file_name,f_len,asctime,dir,DIR_LEN,
     1                      var_in,comment,
     1                      laps_dom_file,DOMAIN_NAME_LEN,imax,jmax,
     1                      n_grids,data,zin,model,grid_spacing,status)

      if (status .eq. -2) GOTO 940
      if (status .eq. -3) GOTO 950
      if (status .eq. -4) GOTO 960
      if (status .eq. -5) GOTO 970
      if (status .eq. -6) GOTO 980
C
C ****  Return normally.
C
      status = ERROR(1)
999   return
C
C ****  Error trapping.
C
930   if (flag .ne. 1)
     1   write(6,*) ' Can only write out ',MAX_GRIDS,' grids.'
      status = ERROR(2)
      goto 999
C
940   if (flag .ne. 1)
     1write(6,*) 'Error opening file to be written to...write aborted.'
      status = ERROR(2)
      goto 999
C
950   if (flag .ne. 1)
     1   write(6,*) 'Error in imax,jmax, or n_grids...write aborted.'
      status = ERROR(2)
      goto 999
C
960   if (flag .ne. 1)
     1   write(6,*) 'Error writing data to file...write aborted.'
      status = ERROR(2)
      goto 999
C
970   if (flag .ne. 1)
     1write(6,*) 'Error writing file header info...write aborted.'
      status = ERROR(2)
      goto 999
C
980   if (flag .ne. 1)
     1   write(6,*) 'Missing one of LAT,LON or AVG grids.'
      status = ERROR(2)
      goto 999
C
      end

