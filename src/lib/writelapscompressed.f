
      subroutine write_laps_compressed(i4time,dir,ext,imax,jmax,
     1   kmax,kdim,var,lvl,lvl_coord,units,comment,data,
     1   istatus)

C**********************************************************************
C
!       implicit  none
C
      integer*4      i4time,               !INPUT I4time of data
     1               i4_valtime,
     1               imax,jmax,kmax,       !INPUT # cols, # rows, # fields
     1               kdim,                 !INPUT K dimension of DATA array
     1               lvl(kdim),            !INPUT Level of each field 
     1               istatus               !OUTPUT

      real*4         data(imax,jmax,kdim)    !INPUT Raw data to be written
      real*4         array(imax*jmax*kdim,2) !LOCAL Compressed array
      character*(*)  dir                     !INPUT Directory to be written to
      character*(*)  ext                     !INPUT File name ext
      character*(*)  var(kdim)               !INPUT 3 letter ID of each field
      character*(*)  lvl_coord(kdim)         !INPUT Vertical coordinate of fields
      character*(*)  units(kdim)             !INPUT units of each field
      character*(*)  comment(kdim)           !INPUT Comments for each field
C
      integer*4      flag,                 !Print flag (1 = off)
     1               i_reftime,            !UNIX time of data
     1               i_valtime,            !UNIX time of data
     1               error(2),
     1               i,j,n7g_nx, n7g_ny,
     1               lgfc,
     1               ldf_len,
     1               fn_length,
     1               var_len,
     1               comm_len,
     1               ext_len,
     1               asc_len,
     1               lvl_coord_len,
     1               units_len,
     1               cdl_path_len,
     1               stat_len,
     1               n_levels,
     1               max_levels,	   !maximum vertical levels
     1               called_from,          !0=FORTRAN, 1=C
     1               append                !0=no, 1=yes
C
      parameter (max_levels=100)
      real*4         pr(max_levels),       !pressures read from get_pres_1d
     1               cdl_levels(max_levels)
C
      character*4    fcst_hh_mm
      character*9    gtime
      character*150  file_name
      character*150  cdl_path
      character*150  static_path
      character*9    laps_dom_file
      character*24   asctime
      character*20   v_g
C
      common         /prt/flag
C
!     include 'lapsparms.cmn'
!     include 'grid_fname.cmn'
C
C-------------------------------------------------------------------------------
C
      error(1)=1
      error(2)=0

C
C ****  Various checks on input data.
C
      if (kmax .gt. kdim) then
        if (flag .ne. 1)
     1write (6,*) 'Illegal K dimension in DATA array...write aborted.'
        istatus=error(2)
        return
      endif
C
      if (imax .ne. n7g_nx) then
        if (flag .ne. 1)
     1write (6,*) 
     1'imax passed in does not match ','laps_dom_file'
     1,'...write aborted.'
        istatus=error(2)
        return
      endif
C
      if (jmax .ne. n7g_ny) then
        if (flag .ne. 1)
     1write (6,*) 
     1'jmax passed in does not match, ',laps_dom_file
     1,'...write aborted.'
        istatus=error(2)
        return
      endif

C ****  Specify file name
C
      call make_fnam_lp(i4time,gtime,istatus)
      if (istatus .ne. 1) then
        write (6,*)
     1'Error converting i4time to file name...write aborted.'
        istatus=error(2)
        return
      endif
C
C **** get actual reftime from gtime...
C
      i_reftime = i4time - 315619200
      i_valtime = i_reftime

C
C ****  Create ascii time variables.
C
      i4_valtime = i_valtime +  315619200
      call cv_i4tim_asc_lp(i4_valtime,asctime,istatus)

      call s_len(ext, ext_len)

      fcst_hh_mm = '0000'

      call cvt_fname_v3(dir,gtime,fcst_hh_mm,ext,ext_len,
     1                  file_name,fn_length,istatus)
      if (istatus .eq. error(2)) goto 930

      called_from = 0    !called from FORTRAN
      append = 0         ! only one analysis time allowed per file

      call s_len(laps_dom_file, ldf_len)
      var_len = len(var(1))
      comm_len = len(comment(1))
      lvl_coord_len = len(lvl_coord(1))
      units_len = len(units(1))
      asc_len = len(asctime)

      ngrids = imax*jmax*kmax
      call runlength_encode(ngrids,n_cmprs_max                ! I
     1                     ,n_cmprs,array,istatus)            ! O

      lun = 65
      call open_lapsprd_file(lun,i4time,ext,istatus)
      if(istatus .ne. 1)goto 940

      write(lun)((array(i,j),j=1,2),i=1,n_cmprs)

      close(lun)
C
C **** write out compressed file
C
      write(6,*) 'laps_dom_file= ',laps_dom_file
      call write_cdf_v3 (file_name,ext,var,comment,asctime,
     1                   cdl_path, static_path, laps_dom_file,
     1                   ldf_len, fn_length,
     1                   ext_len,var_len, 
     1                   comm_len, asc_len, cdl_path_len, 
     1                   stat_len, i_reftime, i_valtime,imax, 
     1                   jmax, kmax, kdim, 
     1                   lvl, data, pr, n_levels, cdl_levels,
     1                   called_from,append, istatus) 
C
      if (istatus .gt. 0) goto 980
      IF (istatus .eq. -2) goto 940
      IF (istatus .eq. -3) goto 950
      IF (istatus .eq. -4) goto 960
      IF (istatus .eq. -5) goto 970
      IF (istatus .eq. -6) goto 990
C
C ****  Return normally.
C
        ISTATUS=ERROR(1)
999     RETURN
C
C ****  Error trapping.
C
920     IF (FLAG .NE. 1) THEN
          write(6,*) 'write_laps ABORTED!'
          write(6,*) ' LAPS will currently only work on a PRESSURE'
     1,' vertical grid'
          if (laps_dom_file(1:lgfc) .eq. 'nest7grid') then
            write(6,*) ' Make sure VERTICAL_GRID is set to PRESSURE'
     1,' in nest7grid.parms'
          else
            write(6,*) ' Make sure VERTICAL_GRID is set to PRESSURE'
     1,' in ',laps_dom_file,'.nl'
          endif
        ENDIF
        ISTATUS=ERROR(2)
        GOTO 999

930     if (flag .ne. 1)
     1    write (6,*) 'file_name variable too short...write aborted.'
        istatus=error(2)
        goto 999
C
940     IF (FLAG .NE. 1)
     1    write (6,*) 'Error opening file to be written to...write abort
     1ed.'
        ISTATUS=ERROR(2)
        GOTO 999
C
950     IF (FLAG .NE. 1)
     1    write (6,*) 'Error in imax,jmax,or n_levels..write aborted'
        ISTATUS=ERROR(2)
        GOTO 999
C
960     IF (FLAG .NE. 1)
     1    write (6,*) 'Error writing data to file...write aborted.'
        ISTATUS=ERROR(2)
        GOTO 999
C
970     IF (FLAG .NE. 1)
     1    write (6,*) 
     1 'Error writing header info into file...write aborted.'
        ISTATUS=ERROR(2)
        GOTO 999
C
980     IF (FLAG .NE. 1)
     1    write (6,*) 
     1 'Some grids not written....could not convert LAPS variables.'
        ISTATUS=ERROR(2)
        GOTO 999
C
990     IF (FLAG .NE. 1)
     1    write (6,*) 
     1 'File already exists for analysis time...write aborted.'
        ISTATUS=ERROR(2)
        GOTO 999
C
        END


        subroutine runlength_encode(ngrids,n_cmprs_max,data   ! I
     1                             ,n_cmprs,array,istatus)    ! O

        real*4 array(n_cmprs_max,2)
        real*4 data(ngrids)

!       Setup for first point
        n_cmprs = 0
        i_count_same = 1

        do i = 2,ngrids-1

            if(data(i) .eq. data(i-1))then    
                i_count_same = i_count_same + 1
            else
                n_cmprs = n_cmprs + 1
                array(n_cmprs,1) = i_count_same
                array(n_cmprs,2) = data(i-1)
                i_count_same = 1
            endif

        enddo ! i

!       Take care of the last point
        i = ngrids

        if(data(i) .eq. data(i-1))then    
            i_count_same = i_count_same + 1
        else
            i_count_same = 1
        endif

        n_cmprs = n_cmprs + 1
        array(n_cmprs,1) = i_count_same
        array(n_cmprs,2) = data(i)

        write(6,*)' End of runlength_encode, number of pts = '
     1           ,n_cmprs,ngrids       
        write(6,*)' Compression ratio = ',float(n_cmprs)/float(ngrids)

        return
        end
