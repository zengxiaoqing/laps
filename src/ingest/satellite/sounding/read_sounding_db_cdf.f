      Subroutine Read_sounder_db_cdf(filename,
     &                               imax,jmax,nch,
     &                               sounding,
     &                               wavelength,
     &                               scalingBias,
     &                               scalingGain,
     &                               northwest_sdr_pixel,
     &                               northwest_sdr_line,
     &                               eastWestCycles,
     &                               eastWestIncs,
     &                               northSouthCycles,
     &                               northSouthIncs,
     &                               frameStartTime,
     &                               lineTimeBegin,lineTimeEnd,
     &                               imc,
     &                               orbitAttitude,
     &                               ndsize_x,ndsize_y,ndsize_ch,
     &                               istatus)
c
c
      implicit none

      integer*4 irec_max
      parameter (irec_max=5000)

      Integer*4  NVARS
      PARAMETER (NVARS=60) !NUMBER OF VARIABLES
C     VARIABLE IDS RUN SEQUENTIALLY FROM 1 TO NVARS= 60
      INTEGER*4 RCODE
C     ****VARIABLES FOR THIS NETCDF FILE****
c
c the following declarations are for testing using the original netCDF read
c routine. Ie., read the entire array all at once
c

      Integer*4   imax,jmax,nch
      Integer*4   sounding(imax,jmax,nch)
      Integer*4   i,j,k,n
      Integer*4   dim_id_y

      REAL*8      wavelength      ( nch )
      Real*4      scalingBias     (jmax,nch)
      Real*4      scalingGain     (jmax,nch)
      Real*4      scaling_rec     (irec_max)

      INTEGER*2   northwest_sdr_pixel
      INTEGER*2   northwest_sdr_line
      INTEGER*4   varid
      Integer*4   ncid
      Integer*4   ncopn
      CHARACTER*1 imc                      (   4)
      REAL*4      imcEnableTime
      INTEGER*2   eastWestCycles                 
      INTEGER*2   eastWestIncs                   
      INTEGER*2   northSouthCycles               
      INTEGER*2   northSouthIncs                 
      REAL*8      frameStartTime                 
      REAL*8      orbitAttitude            (336)
      Real*8      lineTimeBegin(jmax,nch)
      Real*8      lineTimeEnd(jmax,nch)
      Real*8      ltrec(irec_max)

      Integer*4 Istatus
      INTEGER*4 NCVID,ncdid
      Integer*4 ndsize_x(jmax)
      Integer*4 ndsize_y
      Integer*4 ndsize_ch
      Integer*4 ndsize_sb(nch)
      Integer*4 ndsize_ltb(nch)
      Integer*4 ndsize_lte(nch)
      Integer*4 ndsize

      Integer*4 NCNOWRIT
      Integer*4 NVDIM
      Integer*4 NTP,NVS
      Integer*4 LENSTR
      INTEGER*4 START(10)
      INTEGER*4 COUNT(10)
      INTEGER VDIMS(10) !ALLOW UP TO 10 DIMENSIONS
      CHARACTER*31 DUMMY
      CHARACTER*255 filename
C
      istatus = 1

      NCID=NCOPN(filename,
     +NCNOWRIT,RCODE)
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
      call rddata_line_elem(ncid,imax,jmax,nch,
     &ndsize_x,ndsize_y,ndsize_ch,sounding,istatus)
      if(istatus .ne. 1)then
         write(6,*)'Error reading sounding - rddata_img_line_elem'
         return
      endif
c -------------------------------------
C
C    statements to fill wavelength                     
C
      varid = ncvid(ncid,'wavelength',rcode)
      if(rcode.ne.0) return
      CALL NCVINQ(NCID, varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO  20 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
  20  CONTINUE
      CALL NCVGT(NCID, varid,START,COUNT,
     +wavelength                     ,RCODE)
C
C    statements to fill northwest_sdr_pixel            
C
      varid = ncvid(ncid,'northwest_sdr_pixel',rcode)
      if(rcode.ne.0) return
      CALL NCVINQ(NCID, varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO  50 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
  50  CONTINUE
      CALL NCVGT(NCID, varid,START,COUNT,
     +northwest_sdr_pixel            ,RCODE)
C
C    statements to fill northwest_sdr_line             
C
      varid = ncvid(ncid,'northwest_sdr_line',rcode)
      if(rcode.ne.0) return
      CALL NCVINQ(NCID, varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO  60 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
  60  CONTINUE
      CALL NCVGT(NCID, varid,START,COUNT,
     +northwest_sdr_line             ,RCODE)
C
C    statements to fill frameStartTime                 
C
      varid = ncvid(ncid,'frameStartTime',rcode)
      if(rcode.ne.0) return
      CALL NCVINQ(NCID,varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 70 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 70   CONTINUE
      CALL NCVGT(NCID,varid,START,COUNT,
     +frameStartTime                 ,RCODE)
C
C    statements to fill imc                            
C
      varid = ncvid(ncid,'imc',rcode)
      if(rcode.ne.0) return
      CALL NCVINQ(NCID,varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 180 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 180  CONTINUE
      CALL NCVGTC(NCID,varid,START,COUNT,
     +imc                            ,LENSTR,RCODE)
C
C    statements to fill imcEnableTime                  
C
      varid = ncvid(ncid,'imcEnableTime',rcode)
      if(rcode.ne.0) return
      CALL NCVINQ(NCID,varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 260 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 260  CONTINUE
      CALL NCVGT(NCID,varid,START,COUNT,
     +imcEnableTime                  ,RCODE)
C
C    statements to fill eastWestCycles                 
C
      varid = ncvid(ncid,'eastWestCycles',rcode)
      if(rcode.ne.0) return
      CALL NCVINQ(NCID,varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 300 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 300  CONTINUE
      CALL NCVGT(NCID,varid,START,COUNT,
     +eastWestCycles                 ,RCODE)
C
C    statements to fill eastWestIncs                   
C
      varid = ncvid(ncid,'eastWestIncs',rcode)
      if(rcode.ne.0) return
      CALL NCVINQ(NCID,varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 310 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 310  CONTINUE
      CALL NCVGT(NCID,varid,START,COUNT,
     +eastWestIncs                   ,RCODE)
C
C    statements to fill northSouthCycles               
C
      varid = ncvid(ncid,'northSouthCycles',rcode)
      if(rcode.ne.0) return
      CALL NCVINQ(NCID,varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 320 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 320  CONTINUE
      CALL NCVGT(NCID,varid,START,COUNT,
     +northSouthCycles               ,RCODE)
C
C    statements to fill northSouthIncs                 
C
      varid = ncvid(ncid,'northSouthIncs',rcode)
      if(rcode.ne.0) return
      CALL NCVINQ(NCID,varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 330 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 330  CONTINUE
      CALL NCVGT(NCID,varid,START,COUNT,
     +northSouthIncs                 ,RCODE)
C
C    statements to fill orbitAttitude                  
C
      varid = ncvid(ncid,'orbitAttitude',rcode)
      if(rcode.ne.0) return
      CALL NCVINQ(NCID,varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 340 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 340  CONTINUE
      CALL NCVGT(NCID,varid,START,COUNT,
     +orbitAttitude                  ,RCODE)
c
c new: 11-21-96. JRSmart. Retrieve scalinggain and scalingbias. These needed
c to convert counts to useable brightness temps and radiances.
c
      dim_id_y = ncdid(ncid, 'y', rcode)

      varid = ncvid(ncid,'scalingBias',rcode)
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

         call ncdinq(NCID,dim_id_y,dummy,NDSIZE_SB(k),RCODE)
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
         call NCVGT(NCID,varid,START,COUNT,scaling_rec,rcode)
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
      varid = ncvid(ncid,'scalingGain',rcode)
      if(rcode.ne.0)then
         write(6,*)'Error getting scalingGain '
         istatus = -1
      endif

      write(6,*)'Reading ScalingGain'
      do k = 1,nch  !# of channels in sounding database (= 19).

c        write(6,*)'Reading ScalingGain for channel ',k

         call ncdinq(NCID,dim_id_y,dummy,NDSIZE_SB(k),RCODE)
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
         call NCVGT(NCID,varid,START,COUNT,scaling_rec,rcode)
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
      varid = ncvid(ncid,'lineTimeBegin',rcode)
      if(rcode.ne.0)then
         write(6,*)'Error getting lineTimeBegin '
         istatus = -1
      endif
      write(6,*)'Reading lineTimeBegin'
      do k = 1,nch  !# of channels in sounding database (= 19).
         call ncdinq(NCID,dim_id_y,dummy,NDSIZE_LTB(k),RCODE)
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
         call NCVGT(NCID,varid,START,COUNT,ltrec,rcode)
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
      varid = ncvid(ncid,'lineTimeEnd',rcode)
      if(rcode.ne.0)then
         write(6,*)'Error getting lineTimeBegin '
         istatus = -1
      endif
      write(6,*)'Reading lineTimeEnd'
      do k = 1,nch  !# of channels in sounding database (= 19).
         call ncdinq(NCID,dim_id_y,dummy,NDSIZE_LTE(k),RCODE)
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
         call NCVGT(NCID,varid,START,COUNT,ltrec,rcode)
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

      Return
      END
