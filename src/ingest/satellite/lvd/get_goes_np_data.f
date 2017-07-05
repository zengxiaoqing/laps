      subroutine get_goes_np_data
     +                   (i4time_sys,ilaps_cycle_time,NX_L,NY_L
     +                   ,i4time_offset
     +                   ,maxchannel,max_files,nchannels
     +                   ,csat_id,csat_type,chtype
     +                   ,path
     +                   ,image_13,n_ir_elem,n_ir_lines
     +                   ,image_07
     +                   ,image_02,n_vis_elem,n_vis_lines
     +                   ,i4time_goes,istatus)                      ! O

      include 'netcdf.inc'

      character*(*) path
      character*255 filename

      integer nelem, nlines,nf_fid, nf_vid, nf_status
      real image_13(n_ir_elem,n_ir_lines)
      real image_07(n_ir_elem,n_ir_lines)
      real image_02(n_vis_elem,n_vis_lines)

      Character*6  csat_id
      Character*3  csat_type
      Character*3  chtype(maxchannel)
      Character*13 cfname13,cvt_i4time_wfo_fname13

      write(6,*)' Subroutine get_goes_np_data ',n_ir_elem,n_ir_lines

      i4time_goes = ((i4time_sys + 0)
     1            / ilaps_cycle_time) * ilaps_cycle_time + i4time_offset

!     Read IR channel 13 data
      cfname13 = cvt_i4time_wfo_fname13(i4time_goes)
      filename = trim(path)//'/goes_mosaic_'//cfname13//'_13.nc'
      write(6,*)'gnp filename is ',trim(filename)
C
C  Open netcdf File for reading
C
      nf_status=NF_OPEN(filename,NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),filename
        istatus=0
        return
      endif
C
C  Fill all dimension values
C
C
C Get size of nelem
C
      nf_status=NF_INQ_DIMID(nf_fid,'nelem',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim nelem'
      endif
      nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,nelem)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim nelem'
      endif
C
C Get size of nlines
C
      nf_status=NF_INQ_DIMID(nf_fid,'nlines',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim nlines'
      endif
      nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,nlines)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim nlines'
      endif

      call read_goes_np_data(nf_fid,n_ir_elem,n_ir_lines,
     +     csat_id,csat_type,chtype(3),
     +     ilaps_cycle_time, NX_L, NY_L, i4time_earliest,
     +     i4time_latest, lun_out, image_13, istatus)

      nf_status=nf_close(nf_fid)
      
      if(.true.)then

!     Read IR channel 07 data
      cfname13 = cvt_i4time_wfo_fname13(i4time_goes)
      filename = trim(path)//'/goes_mosaic_'//cfname13//'_07.nc'
      write(6,*)'filename is ',trim(filename)
C
C  Open netcdf File for reading
C
      nf_status=NF_OPEN(filename,NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),filename
        istatus=0
        return
      endif
C
C  Fill all dimension values
C
C
C Get size of nelem
C
      nf_status=NF_INQ_DIMID(nf_fid,'nelem',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim nelem'
      endif
      nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,nelem)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim nelem'
      endif
C
C Get size of nlines
C
      nf_status=NF_INQ_DIMID(nf_fid,'nlines',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim nlines'
      endif
      nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,nlines)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim nlines'
      endif

      call read_goes_np_data(nf_fid,n_ir_elem,n_ir_lines,
     +     csat_id,csat_type,chtype(2),
     +     ilaps_cycle_time, NX_L, NY_L, i4time_earliest,
     +     i4time_latest, lun_out, image_07, istatus)

      nf_status=nf_close(nf_fid)

!     Read IR channel 02 data
      cfname13 = cvt_i4time_wfo_fname13(i4time_goes)
      filename = trim(path)//'/goes_mosaic_'//cfname13//'_02.nc'
      write(6,*)'filename is ',trim(filename)
C
C  Open netcdf File for reading
C
      nf_status=NF_OPEN(filename,NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),filename
        istatus=0
        return
      endif
C
C  Fill all dimension values
C
C
C Get size of nelem
C
      nf_status=NF_INQ_DIMID(nf_fid,'nelem',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim nelem'
      endif
      nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,nelem)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim nelem'
      endif
C
C Get size of nlines
C
      nf_status=NF_INQ_DIMID(nf_fid,'nlines',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim nlines'
      endif
      nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,nlines)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim nlines'
      endif

      call read_goes_np_data(nf_fid,n_vis_elem,n_vis_lines,
     +     csat_id,csat_type,chtype(1),
     +     ilaps_cycle_time, NX_L, NY_L, i4time_earliest,
     +     i4time_latest, lun_out, image_02, istatus)

      endif

      nf_status=nf_close(nf_fid)

      write(6,*)' call qc_goes_vis'

      call qc_goes_vis(image_02,n_vis_elem,n_vis_lines)

      istatus = 1

      return
      end
C
C
      subroutine read_goes_np_data(nf_fid, nelem, nlines, 
     +     csat_id,csat_type,chtype,
     +     ilaps_cycle_time, NX_L, NY_L, i4time_earliest,
     +     i4time_latest, lun_out, image, istatus)

      include 'netcdf.inc'
      integer nelem, nlines,nf_fid, nf_vid, nf_status
      real image( nelem, nlines)

!     Declarations for 'write_lvd' call
      integer iwmostanum(recNum)
      character a9time_ob_r(recNum)*9
      logical l_closest_time, l_closest_time_i, l_in_domain
      Character*6  csat_id
      Character*3  csat_type
      Character*3  chtype

      call get_r_missing_data(r_missing_data,istatus)
      if (istatus .ne. 1) then
          write (6,*) 'Error getting r_missing_data'
          return
      endif
       
      image = r_missing_data ! initialize
      
      if(.false.)then
          call read_goes_np_netcdf(nf_fid, nelem, nlines, image)
      else
          nf_status=NF_INQ_VARID(nf_fid,'image',nf_vid)
          if(nf_status.ne.NF_NOERR) then
            print *, NF_STRERROR(nf_status),' for image'
          else
            write(6,*)' call rdblock_line_elem for chtype ',chtype
            call rdblock_line_elem(csat_id,csat_type,chtype,
     &                  nf_fid,nf_vid,nelem,nlines,image,istatus)
          endif
      endif
C
C The netcdf variables are filled - your lvd write call may go here
C
      return
      end
C
C  Subroutine to read the file 
C
      subroutine read_goes_np_netcdf(nf_fid, nelem, nlines, image)
C
      include 'netcdf.inc'
      integer nelem, nlines,nf_fid, nf_vid, nf_status
      integer image( nelem, nlines)




C   Variables of type REAL

C   Variables of type INT
C
C
C     Variable        NETCDF Long Name
C     image         
C
      nf_status=NF_INQ_VARID(nf_fid,'image',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for image'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,image)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for image'
       endif
      endif

C   Variables of type DOUBLE
C


C   Variables of type CHAR
C

      nf_status=nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
      endif

      return
      end

      subroutine qc_goes_vis(image_02,n_vis_elem,n_vis_lines)

      real image_02(n_vis_elem,n_vis_lines)

      iwrite = 0
      nqc = 0

      do ielem = 2,n_vis_lines-1
      do jline = 2,n_vis_lines-1

         ip = ielem+1
         im = ielem-1
         jp = jline+1
         jm = jline-1

         if(image_02(ielem,jline) .eq. 0. .and.
     1      maxval(image_02(im:ip,jm:jp)) .gt. 3000.)then
             nqc = nqc + 1
             if(nqc .le. 100)then
               write(6,*)' VIS QC detected at ',ielem,jline
             endif
             image_02(ielem,jline) = 4095.
         endif         

      enddo ! jline
      enddo ! ielem

      write(6,*)'qc_goes_vis number of QC pixels is ',nqc
      return
      end
