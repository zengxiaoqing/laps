      Subroutine Read_sounder_db_cdf(filename,
     &                               imax,jmax,nch,
     &                               sounding,
     &                               wavelength,
     &                               scalingBias,
     &                               scalingGain,
     &                               northwest_sdr_pixel,
     &                               northwest_sdr_line,
     &                               southeast_sdr_pixel,
     &                               southeast_sdr_line,
     &                               eastWestCycles,
     &                               eastWestIncs,
     &                               northSouthCycles,
     &                               northSouthIncs,
     &                               frameStartTime,
     &                               lineTimeBegin,lineTimeEnd,
     &                               imc,ires_x,ires_y,
     &                               orbitAttitude,
     &                               istatus)
c
c
      implicit none
      include 'netcdf.inc'
      integer irec_max
      parameter (irec_max=5000)

      Integer  NVARS
      PARAMETER (NVARS=60) !NUMBER OF VARIABLES
C     VARIABLE IDS RUN SEQUENTIALLY FROM 1 TO NVARS= 60
      INTEGER RCODE
C     ****VARIABLES FOR THIS NETCDF FILE****
c
c the following declarations are for testing using the original netCDF read
c routine. Ie., read the entire array all at once
c

      Integer   imax,jmax,nch
      Integer   sounding(imax,jmax,nch)
      Integer   i,j,k,n
      Integer   dim_id_y

      REAL*8      wavelength      ( nch )
      Real      scalingBias     (jmax,nch)
      Real      scalingGain     (jmax,nch)
      Real      scaling_rec     (irec_max)

      INTEGER   northwest_sdr_pixel
      INTEGER   northwest_sdr_line
      INTEGER   southeast_sdr_pixel
      INTEGER   southeast_sdr_line
      INTEGER   varid
      CHARACTER*1 imc                      (   4)
      Real      imcEnableTime
      INTEGER   eastWestCycles                 
      INTEGER   eastWestIncs                   
      INTEGER   northSouthCycles               
      INTEGER   northSouthIncs                 
      Integer   ires_x,ires_y
      REAL*8      frameStartTime                 
      REAL*8      orbitAttitude            (336)
      Real*8      lineTimeBegin(jmax,nch)
      Real*8      lineTimeEnd(jmax,nch)
      Real*8      ltrec(irec_max)

      Integer Istatus
      Integer ndsize_sb(nch)
      Integer ndsize_ltb(nch)
      Integer ndsize_lte(nch)
      Integer ndsize, ncid

      Integer NVDIM
      Integer NTP,NVS
      Integer LENSTR
      INTEGER START(10)
      INTEGER COUNT(10)
      INTEGER VDIMS(10) !ALLOW UP TO 10 DIMENSIONS
      CHARACTER*31 DUMMY
      CHARACTER*255 filename
C
      istatus = 1

      RCODE=NF_OPEN(filename,NF_NOWRITE,NCID)
      if(rcode.ne.0)then
         n=index(filename,' ')
         write(6,*)'Error openning netCDF file'
         write(6,*)'filename: ',filename(1:n-1)
         istatus = -1
         return
      endif
c
c code to get dimension size and read individual element of sounding array
c get dimensions for sounding array (x,y,lambda) [lambda is # of wavelengths]
c This code has now been subroutine-ized; rdimg_line_elem_sub.f.
c
c
      call rddata_line_elem(ncid,imax,jmax,nch,sounding,istatus)

      if(istatus .ne. 1)then
         write(6,*)'Error reading sounding - rddata_img_line_elem'
         return
      endif
c -------------------------------------
C
C    statements to fill wavelength                     
C
      rcode=NF_INQ_VARID(ncid,'wavelength',varid)
      if(rcode.ne.0) return
      CALL NCVINQ(NCID, varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO  20 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
  20  CONTINUE
      RCODE=NF_GET_VARA_DOUBLE(NCID,varid,START,COUNT,wavelength)
C
C    statements to fill northwest_sdr_pixel            
C
      rcode=NF_INQ_VARID(ncid,'northwest_sdr_pixel',varid)
      if(rcode.ne.0) return
      CALL NCVINQ(NCID, varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO  50 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
  50  CONTINUE
      RCODE=NF_GET_VARA_INT(NCID,varid,START,COUNT,northwest_sdr_pixel)
C
C    statements to fill northwest_sdr_line             
C
      rcode=NF_INQ_VARID(ncid,'northwest_sdr_line',varid)
      if(rcode.ne.0) return
      CALL NCVINQ(NCID, varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO  60 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
  60  CONTINUE
      RCODE=NF_GET_VARA_INT(NCID,varid,START,COUNT,northwest_sdr_line)
C
C    statements to fill frameStartTime                 
C
      rcode=NF_INQ_VARID(ncid,'frameStartTime',varid)
      if(rcode.ne.0) return
      CALL NCVINQ(NCID,varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 70 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 70   CONTINUE
      RCODE=NF_GET_VARA_DOUBLE(NCID,varid,START,COUNT,frameStartTime)
C
C    statements to fill imc                            
C
      rcode=NF_INQ_VARID(ncid,'imc',varid)
      if(rcode.ne.0) return
      CALL NCVINQ(NCID,varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 180 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 180  CONTINUE
      RCODE= NF_GET_VARA_TEXT(NCID,varid,START,COUNT,imc)
C
C    statements to fill imcEnableTime                  
C
      rcode=NF_INQ_VARID(ncid,'imcEnableTime',varid)
      if(rcode.ne.0) return
      CALL NCVINQ(NCID,varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 260 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 260  CONTINUE
      RCODE=NF_GET_VARA_INT(NCID,varid,START,COUNT,imcEnableTime)
C
C    statements to fill eastWestCycles                 
C
      rcode=NF_INQ_VARID(ncid,'eastWestCycles',varid)
      if(rcode.ne.0) return
      CALL NCVINQ(NCID,varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 300 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 300  CONTINUE
      RCODE=NF_GET_VARA_INT(NCID,varid,START,COUNT,eastWestCycles)
C
C    statements to fill eastWestIncs                   
C
      rcode=NF_INQ_VARID(ncid,'eastWestIncs',varid)
      if(rcode.ne.0) return
      CALL NCVINQ(NCID,varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 310 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 310  CONTINUE
      RCODE=NF_GET_VARA_INT(NCID,varid,START,COUNT,eastWestIncs)
C
C    statements to fill northSouthCycles               
C
      rcode=NF_INQ_VARID(ncid,'northSouthCycles',varid)
      if(rcode.ne.0) return
      CALL NCVINQ(NCID,varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 320 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 320  CONTINUE
      RCODE=NF_GET_VARA_INT(NCID,varid,START,COUNT,northSouthCycles)
C
C    statements to fill northSouthIncs                 
C
      rcode=NF_INQ_VARID(ncid,'northSouthIncs',varid)
      if(rcode.ne.0) return
      CALL NCVINQ(NCID,varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 330 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 330  CONTINUE
      RCODE=NF_GET_VARA_INT(NCID,varid,START,COUNT,northSouthIncs)
C
C    statements to fill orbitAttitude                  
C
      rcode=NF_INQ_VARID(ncid,'orbitAttitude',varid)
      if(rcode.ne.0) return
      CALL NCVINQ(NCID,varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 340 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 340  CONTINUE
      RCODE=NF_GET_VARA_DOUBLE(NCID,varid,START,COUNT,orbitAttitude)
c
c new: 11-21-96. JRSmart. Retrieve scalinggain and scalingbias. These needed
c to convert counts to useable brightness temps and radiances.
c
      dim_id_y = NCDID(ncid, 'y', rcode)

      rcode=NF_INQ_VARID(ncid,'scalingBias',varid)
      if(rcode.ne.0)then
         write(6,*)'Error getting scalingBias '
         istatus = -1
      endif

      if(rcode.ne.0)then
         write(6,*)'Error getting scalingbias id code - returning'
         istatus = -1
         return
      endif

      write(6,*)'Reading ScalingBias'
      do k = 1,nch  !# of channels in sounding database (= 19).

c        write(6,*)'Reading ScalingBias for channel ',k

         call NCDINQ(NCID,dim_id_y,dummy,NDSIZE_SB(k),RCODE)
         if(rcode.ne.0)then
            write(6,*)'Error getting SB dimension - NDSIZE_SB'
         endif

         if(NDSIZE_SB(k).gt.jmax.or.NDSIZE_SB(k).le.0)then
            write(6,*)'SB dimension size > nlines_max or <= 0'
            write(6,*)'ndsize_db(k) = ', ndsize_sb(k)
            istatus = -1
            write(6,*)'Returning to main, istatus = ',istatus
            return
         endif

         COUNT(1)=NDSIZE_SB(k)
         START(2)=k
c
c read record
c
      rcode=NF_GET_VARA_REAL(NCID,varid,START,COUNT,scaling_rec)
         if(rcode.ne.0)then 
            write(6,*)'Error reading scaling database'
            istatus = -1
            write(6,*)'Returning to main, istatus = ',istatus
         endif
c
c load values for all lines into the 3d output array
c
         do i = 1,NDSIZE_SB(k)
            scalingbias(i,k) = scaling_rec(i)
         enddo

      enddo
c
c  statements to fill scalingGain
c
      rcode=NF_INQ_VARID(ncid,'scalingGain',varid)
      if(rcode.ne.0)then
         write(6,*)'Error getting scalingGain '
         istatus = -1
      endif

      write(6,*)'Reading ScalingGain'
      do k = 1,nch  !# of channels in sounding database (= 19).

c        write(6,*)'Reading ScalingGain for channel ',k

         call NCDINQ(NCID,dim_id_y,dummy,NDSIZE_SB(k),RCODE)
         if(rcode.ne.0)then
            write(6,*)'Error getting SG dimension - NDSIZE_SB'
         endif

         if(NDSIZE_SB(k).gt.jmax.or.NDSIZE_SB(k).le.0)then
            write(6,*)'SB dimension size > jmax or <= 0'
            write(6,*)'nsize_sb(k) = ',ndsize_sb(k)
            istatus = -1
            write(6,*)'Returning to main, istatus = ',istatus
            return
         endif

         COUNT(1)=NDSIZE_SB(k)
         START(2)=k
c
c read record
c
      rcode=NF_GET_VARA_REAL(NCID,varid,START,COUNT,scaling_rec)
         if(rcode.ne.0)then
            write(6,*)'Error reading scaling database'
            istatus = -1
            write(6,*)'Returning to main, istatus = ',istatus
         endif
c
c load values for all lines into the 3d output array
c
         do i = 1,NDSIZE_SB(k)
            scalingGain(i,k) = scaling_rec(i)
         enddo

      enddo
c
c 11-14-97 (J.Smart). Acquire lineTimeBegin and lineTimeEnd.
c
      rcode=NF_INQ_VARID(ncid,'lineTimeBegin',varid)
      if(rcode.ne.0)then
         write(6,*)'Error getting lineTimeBegin '
         istatus = -1
      endif
      write(6,*)'Reading lineTimeBegin'
      do k = 1,nch  !# of channels in sounding database (= 19).
         call NCDINQ(NCID,dim_id_y,dummy,NDSIZE_LTB(k),RCODE)
         if(rcode.ne.0)then
            write(6,*)'Error getting LTB dimension - NDSIZE_LTB'
         endif
         if(NDSIZE_LTB(k).gt.jmax.or.NDSIZE_LTB(k).le.0)then
            write(6,*)'LTB dimension size > jmax or <= 0'
            write(6,*)'nsize_ltb(k) = ',ndsize_ltb(k)
            istatus = -1
            write(6,*)'Returning to main, istatus = ',istatus
            return
         endif
         COUNT(1)=NDSIZE_LTB(k)
         START(2)=k
c
c read record, use scaling_rec for the i/o.
c
      rcode=NF_GET_VARA_DOUBLE(NCID,varid,START,COUNT,ltrec)
         if(rcode.ne.0)then
            write(6,*)'Error reading scaling database'
            istatus = -1
            write(6,*)'Returning to main, istatus = ',istatus
         endif
c
c load values for all lines into the 3d output array
c
         do i = 1,NDSIZE_LTB(k)
            lineTimeBegin(i,k) = ltrec(i)
         enddo
      enddo
c -----------
c lineTimeEnd
c
      rcode=NF_INQ_VARID(ncid,'lineTimeEnd',varid)
      if(rcode.ne.0)then
         write(6,*)'Error getting lineTimeBegin '
         istatus = -1
      endif
      write(6,*)'Reading lineTimeEnd'
      do k = 1,nch  !# of channels in sounding database (= 19).
         call NCDINQ(NCID,dim_id_y,dummy,NDSIZE_LTE(k),RCODE)
         if(rcode.ne.0)then
            write(6,*)'Error getting LTE dimension - NDSIZE_LTE'
         endif
         if(NDSIZE_LTE(k).gt.jmax.or.NDSIZE_LTE(k).le.0)then
            write(6,*)'LTE dimension size > jmax or <= 0'
            write(6,*)'nsize_lte(k) = ',ndsize_lte(k)
            istatus = -1
            write(6,*)'Returning to main, istatus = ',istatus
            return
         endif
         COUNT(1)=NDSIZE_LTE(k)
         START(2)=k
c
c read record, use scaling_rec for the i/o.
c
      rcode=NF_GET_VARA_DOUBLE(NCID,varid,START,COUNT,ltrec)
         if(rcode.ne.0)then
            write(6,*)'Error reading scaling database'
            istatus = -1
            write(6,*)'Returning to main, istatus = ',istatus
         endif
c
c load values for all lines into the 3d output array
c
         do i = 1,NDSIZE_LTE(k)
            lineTimeEnd(i,k) = ltrec(i)
         enddo
      enddo
C
C    statements to fill x_resolution and y_resolution
C
      rcode=NF_INQ_VARID(ncid,'x_resolution',varid)
      if(rcode.ne.0) return
      CALL NCVINQ(NCID, varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO  53 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
  53  CONTINUE
      RCODE=NF_GET_VARA_INT(NCID,varid,START,COUNT,ires_x)
C
      rcode=NF_INQ_VARID(ncid,'y_resolution',varid)
      if(rcode.ne.0) return
      CALL NCVINQ(NCID, varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO  54 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
  54  CONTINUE
      RCODE=NF_GET_VARA_INT(NCID,varid,START,COUNT,ires_y)
C
C    statements to fill southeast_sdr_pixel
C
      rcode=NF_INQ_VARID(ncid,'southeast_sdr_pixel',varid)
      if(rcode.ne.0) return
      CALL NCVINQ(NCID, varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO  80 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
  80  CONTINUE
      RCODE=NF_GET_VARA_INT(NCID,varid,START,COUNT,southeast_sdr_pixel)
C
C    statements to fill southeast_sdr_line
C
      rcode=NF_INQ_VARID(ncid,'southeast_sdr_line',varid)
      if(rcode.ne.0) return
      CALL NCVINQ(NCID, varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO  90 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
  90  CONTINUE
      RCODE=NF_GET_VARA_INT(NCID,varid,START,COUNT,southeast_sdr_line)

      Return
      END
c
c get x/y dimensions for satellite data
c
      subroutine get_line_elem_sounder_cdf(filename,nelems,nlines
     &,nchannels,istatus)

      implicit none

      integer nlines
      integer nelems
      integer nchannels
      integer istatus
      integer ln
      integer ncid
      integer dim_id_x
      integer dim_id_y
      integer rcode

      character    filename*(*)
      character*31 dummy

      include 'netcdf.inc'
c
c open file for reading
c
      istatus = 0

      call s_len(filename,ln)
      print*,'open file for x/y/ch dimension read'
      print*,'filename: ',filename(1:ln)
      RCODE=NF_OPEN(filename,NF_NOWRITE,NCID)
      if(rcode.ne.0)then
         print*,'Error openning netCDF file'
         print*,'filename: ',filename(1:ln)
         return
      endif
c
c This is the number of lines
c
      dim_id_y = NCDID(ncid, 'y', rcode)
      if(rcode.ne.0)then
         print*,'Error getting y id code - returning'
         return
      endif
      CALL NCDINQ(NCID, dim_id_y,dummy,nlines,RCODE)
      if(rcode.ne.0)then
         print*,'Error getting y dimension - nlines'
         return
      endif
c
c get x dimension id
c
      dim_id_x = NCDID(ncid, 'x', rcode)
      if(rcode.ne.0)then
         print*,'Error getting x id code - returning'
         return
      endif

      call NCDINQ(NCID,dim_id_x,dummy,nelems,RCODE)
      if(rcode.ne.0)then
         print*,'Error getting x dimension - nelems'
         return
      endif
c
c get wavelength dimension id (this is the number of channels)
c
      dim_id_x = NCDID(ncid, 'wavelength', rcode)
      if(rcode.ne.0)then
         print*,'Error getting wavelength id code - returning'
         return
      endif

      call NCDINQ(NCID,dim_id_x,dummy,nchannels,RCODE)
      if(rcode.ne.0)then
         print*,'Error getting wavelength dimension - nchannels'
         return
      endif

      rcode = nf_close(ncid)
      if(rcode.ne.NF_NOERR) then
        print *, NF_STRERROR(rcode)
        print *,'nf_close'
        return
      endif

      istatus = 1

      return
      end
