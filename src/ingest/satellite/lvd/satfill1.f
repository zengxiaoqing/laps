      subroutine satfill1(csat_id,csat_type,
     &               i4time_data_in,smsng,
     &               c_type,
     &               n_lines,n_elems,
     &               cdirpath,
     &               r_missing_data,
     &               image_main,
     &               percent_missing,
     &               istatus)
c
c
c
      implicit none
c
      include 'netcdf.inc'

      integer   n_lines
      integer   n_elems
      integer   i,n,j,nt
      integer   record

      REAL      la1
      REAL      lo1
      REAL      dx,dy
      REAL      latin
      REAL      lov
      INTEGER   validTime
      Integer   nx2,ny2
      integer   rcode,ncid

      character*13  cfname13
      character*13  cvt_i4time_wfo_fname13
      character*9   cfname
      character*3   csat_type
      character*3   c_type
      character*6   csat_id
      character*200 cdirpath
      character*200 c_filename
      character*255 pathname

      real    image_main(n_elems,n_lines)
      real    image_prev(n_elems,n_lines)
      real    latitude(n_elems,n_lines)
      real    longitude(n_elems,n_lines)
      real    r_missing_data
      real    percent_missing
      real    smsng
      real    scale_img

      integer bcnt
      integer itot
      integer istatus
      integer lstatus
      integer io_status
      integer iqstatus
      integer i4time_data_in
      integer i4time_needed_in
      integer i4time_nearest
c
c -------------------------------------------------------
c
      istatus = 1

      nt=index(c_type,' ')-1
      if(nt.le.0)nt=3

      bcnt=0
      itot=n_lines*n_elems
      do j = 1,n_lines
         do i = 1,n_elems
            if(image_main(i,j).eq.r_missing_data)then
               bcnt=bcnt+1
            endif
         enddo
      enddo
 
      percent_missing=float(bcnt)/float(itot)
c
c if any bad/missing data points found then acquire image from
c 15 minutes previous
c
      write(6,*)'    Num of msng pix in satfill (current image): ',bcnt
      write(6,*)'    percent missing: ',percent_missing
c
c only try to fill current image with pixels from previous (no more than 15 min)
c image if less than 25% of current image is bad.
c
      if(percent_missing.gt.0.0.and.percent_missing.le.0.25)then

         i4time_needed_in=i4time_data_in-900
         n=index(cdirpath,' ')

         if(csat_type.eq.'wfo')then
            cfname13 = cvt_i4time_wfo_fname13(i4time_needed_in)
            pathname=cdirpath(1:n-1)//cfname13(1:11)//'*'
         else
            call make_fnam_lp(i4time_needed_in,cfname,lstatus)
          pathname=cdirpath(1:n-1)//cfname(1:8)//'*_'//c_type(1:nt)
         endif

         call get_file_time(pathname 
     1          ,i4time_needed_in,i4time_nearest)

         if(abs(i4time_needed_in-i4time_nearest).le.181)then
            write(*,*)'    Found file for qc fill'
            if(csat_type.eq.'wfo')then
               cfname13=cvt_i4time_wfo_fname13(i4time_nearest)
               c_filename=cdirpath(1:n-1)//cfname13
            else
               call make_fnam_lp(i4time_nearest,cfname,lstatus)
               c_filename=cdirpath(1:n-1)//cfname//'_'//c_type(1:nt)
               rcode=NF_OPEN(c_filename,NF_NOWRITE,NCID)
               if(rcode.ne.nf_noerr) return
               call get_cdf_dims(ncid,record,nx2,ny2,istatus)
               record=1
               if(istatus.eq.1)then
                  print*,'Error reading cdf dimensions'
                  return
               endif
 
            endif
            n=index(c_filename,' ')
            write(6,*)'    Filename: ',c_filename(1:n-1)

            call readcdf(csat_id,csat_type,c_type,
!    &                   nx2,ny2,
     &                   record,
     &                   n_elems,n_lines,
     &                   image_prev,scale_img,
     &                   latitude,
     &                   longitude,
     &                   la1,lo1,
     &                   Dx,Dy,
     &                   Latin,Lov,
     &                   validTime,
     &                   ncid,
     &                   io_status)
            if(io_status .ne. 1)then
               write(6,*)'    Error getting image for qc'
               write(6,*)'    Filename: ',c_filename(1:n-1)
               goto 998
            endif
c
c The missing data flag is r_missing_data!
c First examine the previous image for missing data.
c Get missing satellite data value.
c
            call  set_missing_sat(csat_id,csat_type,c_type,
     &               image_prev,n_elems,n_lines,
     &               smsng,r_missing_data,
     &               istatus)
            if(istatus .lt. 0)then
               write(6,*)'    Found missing data in set_missing_sat'
               write(6,*)'    Num of msng pix in prev img',abs(istatus)
               goto 998
            else
               write(6,*)'    No missing pixels in previous image'
               write(6,*)'    using previous image in image_compare'
            endif
c
c use following subroutine to compare the main and previous images
c for qcfill.
c
            call image_compare(n_lines,n_elems,r_missing_data,
     &                         image_main,image_prev,
     &                         iqstatus)

            if(iqstatus.gt.0)then
               write(6,*)'    Sat pixels filled in image_compare'
               write(6,*)'    = ',iqstatus
               write(6,*)'    recompute percent_missing'
c
c recompute the percent missing
c
               bcnt=0
               itot=0
               do j = 1,n_lines
               do i = 1,n_elems
                  itot=itot + 1
                  if(image_main(i,j).eq.r_missing_data)then
                     bcnt=bcnt+1
                  endif
               enddo
               enddo

               percent_missing=float(bcnt)/float(itot)

            elseif(iqstatus.lt.0)then
               write(6,*)'    Error filling points in image_compare'
               goto 1000
            else
               write(6,*)'    No pixels filled in image_compare'
            endif

         else

            write(6,*)'    No previous image exist for fill routine'

         endif

      elseif(percent_missing.ge.0.25)then

         write(6,*)'    More than 25% of current image missing!'
         write(6,*)'    ***************************************'

      else

         write(6,*)'    No filling required for image_main'

      endif

      goto 1000

998   write(6,*)'Not using prev image to fill'

1000  return
      end
