
       subroutine read_laps_compressed(i4time,dir,ext,
     1                    iimax,jjmax,kdim,
     1                    var_req,lvl_req,lvl_coord_req,
     1                    units_req_o,comment_req_o,data_o,istatus)

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

      integer nf
      parameter (nf=3)

      integer       i4time,               !INPUT I4time of data
     1                iimax,jjmax,          !INPUT # cols, # rows
     1                kdim,                 !INPUT K dimension of DATA array
     1                lvl_req(kdim*nf),     !INPUT Requested levels
     1                istatus               !OUTPUT
      character*(*)   dir                   !INPUT Directory to read data from
      character*(*)   ext                   !INPUT File name ext 
      character*(*)   var_req(kdim)         !INPUT 3 letter ID of requested fields
      character*(*)   lvl_coord_req(kdim)   !OUTPUT Vertical coordinate of fields

      character*(*)   units_req_o(kdim)     !OUTPUT Units of requested fields
      character*(10)  units_req_l(kdim,nf)  !LOCAL  Units of requested fields

      character*(*)   comment_req_o(kdim)   !OUTPUT Comments for requested fields
      character*(125) comment_req_l(kdim,nf)!LOCAL  Comments for requested fields

      real        data_o(iimax,jjmax,kdim)               !OUTPUT data
!     real        data_l(iimax,jjmax,kdim,nf)            !LOCAL data

      real,    allocatable, dimension(:,:,:,:) :: data_l !LOCAL data

      integer, allocatable, dimension(:) :: array1       !LOCAL Compressed array
      real,    allocatable, dimension(:) :: array2       !LOCAL Compressed array
C
      integer fn_length,
     1          flag,                   !Print flag (1 = off)
     1          error(3),
     1          comm_len,
     1          ext_len
  
C
      character*5       fcst_hh_mm
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

C **** get actual reftime from gtime...
C
      comm_len = len(comment_req_o(1))
C
C **** read in compressed file
C
      lun = 65

      itry = 0 ! retry if file is only partially written out
5     if(.false.)then ! old way
          call open_lapsprd_file_read(lun,i4time,ext,istatus)
          if(istatus .ne. 1)goto 950
      else            ! new way using 'dir'

C fcst_hh_mm: Hard wired as a place holder - will be used in filename only
C   if read_laps_compressed is called on lga, lgb, fua, fsf, ram, rsf
C To fix this, call read_laps instead
          fcst_hh_mm = '0000'

          call cvt_fname_v3(dir,gtime,fcst_hh_mm,ext,ext_len,
     1                      file_name,fn_length,istatus)
          if (istatus .eq. error(2)) goto 930

          open(lun,file=file_name(1:fn_length),status='old',err=950)
      endif

10    read(lun,*,err=11,end=11)kkdim
      go to 12
11    write(6,*)' Warning: could not read kkdim'
      call sleep(5)
      itry = itry + 1
      if(itry .le. 5)then
          close(lun) 
          write(6,*)' Retrying the read'
          go to 5 
      else
          go to 960
      endif

12    if(kkdim .ne. kdim * nf)then
          write(6,*)kkdim,kdim*nf
          go to 980
      endif

      write(6,*)' Reading comments, length = ',comm_len
      do l = 1,nf
          do k = 1,kdim
              read(lun,1)comment_req_l(k,l)
              ! Temporarily comment this following write statement out for saving output time and file July 2013
              ! Steve and Yuanfu will add a parameter later controling the output.
              ! write(6,1)comment_req_l(k,l)
1             format(a)
          enddo ! k
      enddo ! l

      read(lun,*)n_cmprs

      allocate(array1(n_cmprs),STAT=istat_alloc)
      if(istat_alloc .ne. 0)then
          write(6,*)' ERROR: Could not allocate array1'
          stop
      endif

      allocate(array2(n_cmprs),STAT=istat_alloc)
      if(istat_alloc .ne. 0)then
          write(6,*)' ERROR: Could not allocate array2'
          stop
      endif

      icheck_sum = 0

      do i = 1,n_cmprs
20        read(lun,*,err=21,end=21)array1(i),array2(i)
          go to 22
21        write(6,*)' Warning: could not read array1/array2'
          call sleep(5)
          itry = itry + 1
          if(itry .le. 5)then
              close(lun)
              write(6,*)' Retrying the read'
              go to 5 
          else
              go to 960
          endif

22        icheck_sum = icheck_sum + array1(i)
      enddo ! i

      close(lun)

      ngrids = iimax*jjmax*kdim*nf

      if(icheck_sum .ne. ngrids)then
          write(6,*)icheck_sum, ngrids
          go to 980
      endif

!     Decode the data 
      write(6,*)' Decoding the data'

      n_cmprs_max = ngrids

      allocate(data_l(iimax,jjmax,kdim,nf),STAT=istat_alloc)
      if(istat_alloc .ne. 0)then
          write(6,*)' ERROR: Could not allocate data_l'
          stop
      endif

      call runlength_decode(ngrids,n_cmprs,array1,array2        ! I
     1                     ,data_l                              ! O
     1                     ,istatus)                            ! O
      deallocate(array1)
      deallocate(array2)
      if(istatus .ne. 1)goto 970

!     Position the data into the proper 3-D array
      if(var_req(1) .eq. 'REF')then
          ifield = 1
      elseif(var_req(1) .eq. 'VEL')then
          ifield = 2
      elseif(var_req(1) .eq. 'NYQ')then
          ifield = 3
      else
          write(6,*)' ERROR: unknown variable requested ',var_req
          goto930
      endif

      do k = 1,kdim
!         Assign Comment
          comment_req_o(k) = comment_req_l(k,ifield)

!         Assign Level Coord?

!         Assign Units?

!         do i = 1,iimax
!         do j = 1,jjmax
!             data_o(i,j,k) = data_l(i,j,k,ifield)
!         enddo ! j
!         enddo ! i

      enddo ! k      

      data_o(:,:,:) = data_l(:,:,:,ifield)
C
C ****  Return normally.
C
        istatus=error(1)
999     if(allocated(data_l))deallocate(data_l)
        return
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
960     if (flag .ne. 1)
     1    write (6,*) 'Error reading compressed file...read aborted.'       
        istatus=error(2)
        goto 999
C
970     if (flag .ne. 1)
     1    write (6,*) 'Error decoding data...read aborted.'
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


        subroutine runlength_decode(ngrids,n_cmprs,array1,array2        ! I
     1                     ,data                                        ! O
     1                     ,istatus)                                    ! O

        integer array1(n_cmprs)
        real array2(n_cmprs)
        real data(ngrids)

!       Setup for first point
        i_end = 0

        do i = 1,n_cmprs ! Loop through array list
            i_start = i_end + 1
            i_end = i_start + array1(i) - 1

            do ii = i_start,i_end
                data(ii) = array2(i)
            enddo ! ii

        enddo ! i

        if(i_end .ne. ngrids)then
            write(6,*)' Error in runlength_decode',i_end,ngrids
            istatus = 0
            return
        endif

        istatus = 1
        return
        end
