      subroutine wrt_laps_static (dir, laps_dom_file, imax, jmax, 
     1                            n_grids, dx, dy, lov, latin1, 
     1                            latin2, origin, var, comment, 
     1                            data, model, grid_spacing,
     1                            map_proj, status)

      implicit none

C This routine calls functions in static_routines.c.  To create the new
C static file, a cdl file named <laps_dom_file>.cdl (ie. nest7grid.cdl)
C must be located in the LAPS cdl directory and is written out to the
C directory specified in "dir".

C Passed in variables
      character      dir*(*),laps_dom_file*(*)
      integer*4      imax, jmax, n_grids
      real*4         dx, dy, lov, latin1, latin2
      character*(*)  origin
      character*(*)  var(n_grids)
      character*(*)  comment(n_grids)
      real*4         data(imax,jmax,n_grids)
      character*(*)  model
      real*4         grid_spacing
      character*(*)  map_proj
      integer*4      status

C Local variables

      integer*4      i4time_now_gg, dom_len, map_len
      integer*4      origin_len, unixtime, nx_lp, ny_lp
      integer*4      var_len, com_len, cdl_dir_len, model_len 
      character*132  file_name, cdl_dir
      character*24   asctime
      character*30   map_projection
      
      integer*4      i4time,
     1               ERROR(2),
     1               flag

      integer*4      f_len, len_asc

      COMMON         /PRT/flag

      ERROR(1)=1
      ERROR(2)=0

C  BEGIN SUBROUTINE

      call make_static_fname(dir,laps_dom_file,file_name,f_len,status)

      i4time = i4time_now_gg()
      call cv_i4tim_asc_lp(i4time,asctime,status)
      unixtime = i4time - 315619200

c get NX_L and NY_L from nest7grid.parms
      call get_grid_dim_xy(nx_lp, ny_lp, status)

      call s_len(laps_dom_file, dom_len)
      var_len = len(var(1))
      com_len = len(comment(1))
      origin_len = len(origin)
      model_len = len(model)
      len_asc = len(asctime)
      map_len = len(map_projection)
      map_projection = ' '
      if (map_proj .eq. 'plrstr') 
     1  map_projection = 'polar stereographic'
      if (map_proj .eq. 'lambrt' ) then
        if (latin1 .eq. latin2) then
          map_projection = 'tangential lambert conformal'
        else
          map_projection = 'secant lambert conformal'
        endif
      endif
      if (map_proj .eq. 'merctr') map_projection = 'mercator'

      if (lov .lt. 0.0) lov = 360.0 + lov
     
      call get_directory('cdl',cdl_dir, cdl_dir_len)
      
      print*,'map_len ', map_len

      file_name(132:132) = char(0)
      cdl_dir(132:132) = char(0)
      asctime(24:24) = char(0)
      map_projection(30:30) = char(0)


      call write_cdf_static(file_name,f_len,asctime,cdl_dir
     1     ,cdl_dir_len,var,var_len,comment,com_len
     1     ,laps_dom_file,dom_len,imax,jmax,n_grids
     1     ,nx_lp,ny_lp,data,model,grid_spacing,dx
     1     ,dy,lov,latin1,latin2,origin,origin_len
     1     ,map_projection,map_len,unixtime, status)

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
     1write(6,*) 'x and y values in nest7grid.cdl do not match imax
     1and jmax in nest7grid.parms...write aborted.'
      status = ERROR(2)
      goto 999
C
      end

