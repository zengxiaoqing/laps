
      subroutine write_laps_compressed(i4time,dir,ext,imax,jmax,
     1   kmax,kdim,var,lvl,lvl_coord,units,comment,data,
     1   istatus)

C**********************************************************************
C
!       implicit  none
C
      integer      i4time,               !INPUT I4time of data
     1               i4_valtime,
     1               imax,jmax,kmax,       !INPUT # cols, # rows, # fields
     1               kdim,                 !INPUT K dimension of DATA array
     1               lvl(kdim),            !INPUT Level of each field 
     1               istatus               !OUTPUT

      real         data(imax,jmax,kdim)    !INPUT Raw data to be written

      integer, allocatable, dimension(:) :: array1 !LOCAL Compressed array
      real,    allocatable, dimension(:) :: array2 !LOCAL Compressed array

      character*(*)  dir                     !INPUT Directory to be written to
      character*(*)  ext                     !INPUT File name ext
      character*(*)  var(kdim)               !INPUT 3 letter ID of each field
      character*(*)  lvl_coord(kdim)         !INPUT Vertical coordinate of fields
      character*(*)  units(kdim)             !INPUT units of each field
      character*(*)  comment(kdim)           !INPUT Comments for each field
C
      integer      flag,                 !Print flag (1 = off)
     1               i_reftime,            !UNIX time of data
     1               i_valtime,            !UNIX time of data
     1               error(2),
     1               i,j,n7g_nx, n7g_ny,
     1               lgfc,
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
      real         pr(max_levels),       !pressures read from get_pres_1d
     1               cdl_levels(max_levels)
c
      logical l_check_encoding
C
      character*5    fcst_hh_mm
      character*9    gtime
      character*150  file_name
      character*150  cdl_path
      character*150  static_path
      character*24   asctime
      character*20   v_g
C
      common         /prt/flag
C
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

C fcst_hh_mm: Hard wired as a place holder - will be used in filename  only
C   if write_laps_compressed is called on lga, lgb, fua, fsf, ram, rsf
C To fix this, call write_laps instead
      fcst_hh_mm = '0000'

      call cvt_fname_v3(dir,gtime,fcst_hh_mm,ext,ext_len,
     1                  file_name,fn_length,istatus)
      if (istatus .eq. error(2)) goto 930

      called_from = 0    !called from FORTRAN
      append = 0         ! only one analysis time allowed per file

      var_len = len(var(1))
      comm_len = len(comment(1))
      lvl_coord_len = len(lvl_coord(1))
      units_len = len(units(1))
      asc_len = len(asctime)

!     Allocate arrays then do run-length encoding
!     n_cmprs_max = imax*jmax*kdim
      n_cmprs_max = 2000000  
 800  write(6,*)' Allocate arrays with size of ',n_cmprs_max

      allocate( array1(n_cmprs_max), STAT=istat_alloc )
      if(istat_alloc .ne. 0)then
          write(6,*)' ERROR: Could not allocate array1'
          write(6,*)' Try reducing n_cmprs_max from ',n_cmprs_max
          goto 950
      endif

      allocate( array2(n_cmprs_max), STAT=istat_alloc )
      if(istat_alloc .ne. 0)then
          write(6,*)' ERROR: Could not allocate array2'
          write(6,*)' Try reducing n_cmprs_max from ',n_cmprs_max
          goto 950
      endif

      ngrids = imax*jmax*kdim

      call runlength_encode(ngrids,n_cmprs_max,data           ! I
     1                     ,n_cmprs,array1,array2,istatus)    ! O
      if(istatus .eq. -1)then ! try increasing array allocations
          deallocate(array1)
          deallocate(array2)
          n_cmprs_max = n_cmprs_max * 2
          goto 800
      elseif(istatus .ne. 1)then
          goto 950
      endif

      l_check_encoding = .false.
      if(l_check_encoding)then ! Just for debugging purposes
          call runlength_decode(ngrids,n_cmprs,array1,array2        ! I
     1                         ,data                                ! O
     1                         ,istatus)                            ! O
          if(istatus .ne. 1)then
              write(6,*)
     1            ' 1st Decoding test of compressed data unsuccessful'
              goto 980
          else
              write(6,*)
     1            ' 1st Decoding test of compressed data successful'
          endif
      endif
C
C **** write out compressed file
C
      lun = 65
      call open_lapsprd_file(lun,i4time,ext,istatus)
      if(istatus .ne. 1)goto 940

      write(lun,*)kdim

      do k = 1,kdim
          write(lun,1)comment(k)
1         format(a)
      enddo ! k

      write(lun,*)n_cmprs

      icheck_sum = 0

      do i = 1,n_cmprs
          write(lun,*)array1(i),array2(i)
          icheck_sum = icheck_sum + array1(i)
      enddo ! i

      close(lun)

      deallocate(array1)
      deallocate(array2)

!     Second "internal" checksum test
      if(icheck_sum .ne. ngrids)then
          write(6,*)' 2nd checksum test discrepancy: '
     1             ,icheck_sum, ngrids
          go to 980
      endif
C
C ****  Return normally.
C
        ISTATUS=ERROR(1)
999     RETURN
C
C ****  Error trapping.
C
920     IF (FLAG .NE. 1) THEN
          write(6,*) ' write_laps_compressed ABORTED!'
          write(6,*) ' LAPS will currently only work on a PRESSURE'
     1              ,' vertical grid'
          write(6,*) ' Make sure VERTICAL_GRID is set to PRESSURE'
     1              ,' in nest7grid.parms'
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
     1    write (6,*) 'Error in/near runlength_encode...write aborted'       
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
     1   'Checksum error with runlength encoding'
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


        subroutine runlength_encode(ngrids,n_cmprs_max,data           ! I
     1                             ,n_cmprs,array1,array2,istatus)    ! O

        integer array1(n_cmprs_max)
        real array2(n_cmprs_max)
        real data(ngrids)

!       Setup for first point
        n_cmprs = 0
        i_count_same = 1

        do i = 2,ngrids-1

            if(data(i) .eq. data(i-1))then    
                i_count_same = i_count_same + 1
            else
                n_cmprs = n_cmprs + 1
                if(n_cmprs .le. n_cmprs_max)then
                    array1(n_cmprs) = i_count_same
                    array2(n_cmprs) = data(i-1)
                    i_count_same = 1
                else
                    write(6,*)' ERROR, increase n_cmprs_max',n_cmprs_max
                    istatus = -1
                    return
                endif
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
        if(n_cmprs .le. n_cmprs_max)then
            array1(n_cmprs) = i_count_same
            array2(n_cmprs) = data(i)
        else
            write(6,*)' ERROR, increase n_cmprs_max',n_cmprs_max
            istatus = -1
            return
        endif

        write(6,*)' End of runlength_encode, number of pts = '
     1           ,n_cmprs,ngrids       
        write(6,*)' Compression ratio = ',float(n_cmprs)/float(ngrids)

        istatus = 1
        return
        end
