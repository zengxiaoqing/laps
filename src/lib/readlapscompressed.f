
       subroutine read_laps_compressed(i4time,dir,ext,
     1                    iimax,jjmax,kkmax,kdim,
     1                    var_req,lvl_req,lvl_coord_req,
     1                    units_req,comment_req,data,istatus)

C**********************************************************************
C
C      Subroutine READ_LAPS_COMPRESSED
C
C      Author:    Steve Albers
C
C      Reads data requested by arrays VAR_REQ and LVL_REQ for the
C      I4TIME, DIR and EXT specified.  Returns LVL_COORD-REQ,
C      UNITS_REQ, COMMENT_REQ, DATA AND ISTATUS
C
C      Based on 'READ_LAPS_DATA' by Linda Wharton
C
C**********************************************************************
C
C
      integer*4       i4time,              !INPUT I4time of data
     1                iimax,jjmax,kkmax,   !INPUT # cols, # rows, # fields
     1                kdim,                !INPUT K dimension of DATA array
     1                lvl_req(kdim),       !INPUT Requested levels
     1                istatus              !OUTPUT
      character*(*)   dir                  !INPUT Directory to read data from
      character*(*)   ext                  !INPUT File name ext 
      character*(*)   var_req(kdim)        !INPUT 3 letter ID of requested fields
      character*(*)   lvl_coord_req(kdim)  !OUTPUT Vertical coordinate of fields
      character*(*)   units_req(kdim)      !OUTPUT Units of requested fields
      character*(*)   comment_req(kdim)    !OUTPUT Comments for requested fields
      real*4        data(iimax,jjmax,kdim) !OUTPUT data
      real*4        array(iimax*jjmax*kdim,2) !LOCAL Compressed array
C
      integer*4 fn_length,
     1          i_reftime,              !UNIX time of data
     1          i_valtime,              !UNIX time of data
     1          flag,                   !Print flag (1 = off)
     1          error(3),
     1          called_from,
     1          var_len,
     1          comm_len,
     1          ext_len,
     1          lvl_coord_len,
     1          units_len
  
C
      character*4       fcst_hh_mm
      character*9       gtime
      character*150     file_name
C
      common            /prt/flag
C
C-------------------------------------------------------------------------------
C
      error(1)=1
      error(2)=0
      error(3)=-2
C
C ****  Create file name.
C
      call make_fnam_lp(i4time,gtime,istatus)
      if (istatus .ne. 1) then
        write (6,*) 
     1'Error converting i4time to file name...read aborted.'
        istatus=error(2)
        return
      endif

      call s_len(ext, ext_len)

      fcst_hh_mm = '0000'

      call cvt_fname_v3(dir,gtime,fcst_hh_mm,ext,ext_len,
     1                  file_name,fn_length,istatus)
      if (istatus .eq. error(2)) goto 930
C
C **** get actual reftime from gtime...
C
      i_reftime = i4time - 315619200
      i_valtime = i_reftime

      called_from = 0   !called from FORTRAN

      var_len = len(var_req(1))
      comm_len = len(comment_req(1))
      lvl_coord_len = len(lvl_coord_req(1))
      units_len = len(units_req(1))
C
C **** read in compressed file
C
      write(6,*) 'laps_dom_file= ',laps_dom_file

      lun = 65
      call open_lapsprd_file(lun,i4time,ext,istatus)
      if(istatus .ne. 1)goto 950

      read(lun,*)kdim

      do k = 1,kdim
          read(lun,*)comment_req(k)
      enddo ! k

      read(lun,*)n_cmprs

      icheck_sum = 0

!     read(lun,*)((array(i,j),j=1,2),i=1,n_cmprs)
      do i = 1,n_cmprs
          read(lun,*)array(i,1),array(i,2)
          icheck_sum = icheck_sum + nint(array(i,1))
      enddo ! i

      close(lun)

      if(icheck_sum .ne. ngrids)then
          write(6,*)icheck_sum, ngrids
          go to 980
      endif
C
C ****  Return normally.
C
        istatus=error(1)
999     return
C
C ****  Error trapping.
C
930     if (flag .ne. 1)
     1    write (6,*) 'file_name variable too short...read aborted.'
        istatus=error(2)
        goto 999
C
950     if (flag .ne. 1)
     1    write (6,*) 'Error opening compressed file...read aborted.'       
        istatus=error(2)
        goto 999
C
970     if (flag .ne. 1)
     1    write (6,*) 'Error retrieving data...read aborted.'
        istatus=error(2)
        goto 999
C
980     if (flag .ne. 1)
     1    write (6,*) 'Error in array dimensioning...read aborted.'
        istatus=error(2)
        goto 999
C
990     if (flag .ne. 1)
     1    write (6,*) 'Error in version, not a valid LAPS file... '
     1                ,'read aborted.'
        istatus=error(2)
        goto 999
C
992     continue
        istatus=error(2)
        goto 999
C
        END


        subroutine runlength_decode(ngrids,n_cmprs_max,data   ! I
     1                             ,n_cmprs,array,istatus)    ! O

!       Still needs to be reworked from original 'encode' routine

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
