      subroutine rd_gvarimg_cdf_header(filename,
     &                              northwest_vis_pixel,
     &                              northwest_vis_line,
     &                              southeast_vis_pixel,
     &                              southeast_vis_line,
     &                              eastWestCycles,
     &                              eastWestIncs,
     &                              northSouthCycles,
     &                              northSouthIncs,
     &                              x_resolution,
     &                              y_resolution,
     &                              frameStartTime,
     &                              imc,
     &                              orbitAttitude,
     &                              Nx4,Ny4,
     &                              x_step,y_step,
     &                              istatus)
c
      PARAMETER (NVARS=67) !NUMBER OF VARIABLES
      include 'netcdf.inc'
c
      INTEGER RCODE
      INTEGER VARID

C     ****VARIABLES FOR THIS NETCDF FILE****
C
      CHARACTER*1 wavelength                     (   1,   3)
      INTEGER   northwest_vis_pixel            (   1)
      INTEGER   northwest_vis_line             (   1)
      INTEGER   southeast_vis_pixel            (   1)
      INTEGER   southeast_vis_line             (   1)
      INTEGER   x_step,y_step
      REAL*8      frameStartTime                 
      CHARACTER*1 elem_dim                       (   1,   8)
      CHARACTER*1 line_dim                       (   1,   8)
      CHARACTER*1 imc                            (   4)
      INTEGER   eastWestCycles                 
      INTEGER   eastWestIncs                   
      INTEGER   northSouthCycles               
      INTEGER   northSouthIncs                 
      REAL*8      orbitAttitude                  (   1, 336)
c     CHARACTER*1 x_dim                          (   1,   8)
c     CHARACTER*1 y_dim                          (   1,   8)
      INTEGER   Nx                             (   1)
      INTEGER   Ny                             (   1)
      REAL        Lap                            (   1)
      REAL        Lop                            (   1)
      Integer   x_resolution
      Integer   y_resolution
      REAL        dx                             (   1)
      REAL        dy                             (   1)
      character   filename*255
c
c*************************************
c
      INTEGER START(10)
      INTEGER COUNT(10)
      INTEGER Nx4,Ny4
      INTEGER VDIMS(10) !ALLOW UP TO 10 DIMENSIONS
      CHARACTER*31 DUMMY
c
c ======  START =======
c
      RCODE=NF_OPEN(filename,NF_NOWRITE,NCID)

C
C    statements to fill wavelength                     
C
       rcode=NF_INQ_VARID(ncid,'wavelength',varid)
      if(rcode.ne.0)return
      CALL NCVINQ(NCID, varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO  90 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
  90  CONTINUE
      RCODE= NF_GET_VARA_TEXT(NCID, varid,START,COUNT,wavelength)
C
C    statements to fill northwest_vis_pixel            
C
       rcode=NF_INQ_VARID(ncid,'northwest_vis_pixel',varid)
      if(rcode.ne.0)return
      CALL NCVINQ(NCID,varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 100 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 100  CONTINUE
      RCODE=NF_GET_VARA_INT(NCID,varid,START,COUNT,northwest_vis_pixel)
C
C    statements to fill northwest_vis_line             
C
       rcode=NF_INQ_VARID(ncid,'northwest_vis_line',varid)
      if(rcode.ne.0)return
      CALL NCVINQ(NCID,varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 110 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 110  CONTINUE
      RCODE=NF_GET_VARA_INT(NCID,varid,START,COUNT,northwest_vis_line)
C
C    statements to fill southeast_vis_pixel            
C
       rcode=NF_INQ_VARID(ncid,'southeast_vis_pixel',varid)
      if(rcode.ne.0)return
      CALL NCVINQ(NCID,varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 120 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 120  CONTINUE
      RCODE=NF_GET_VARA_INT(NCID,varid,START,COUNT,southeast_vis_pixel)
C
C    statements to fill southeast_vis_line             
C
       rcode=NF_INQ_VARID(ncid,'southeast_vis_line',varid)
      if(rcode.ne.0)return

      CALL NCVINQ(NCID,varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 130 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 130  CONTINUE
      RCODE=NF_GET_VARA_INT(NCID,varid,START,COUNT,southeast_vis_line)
C
C    statements to fill frameStartTime                 
C
       rcode=NF_INQ_VARID(ncid,'frameStartTime',varid)
      if(rcode.ne.0)return

      CALL NCVINQ(NCID,varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 150 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 150  CONTINUE
      RCODE=NF_GET_VARA_DOUBLE(NCID,varid,START,COUNT,frameStartTime)
C
C    statements to fill elem_dim                       
C
       rcode=NF_INQ_VARID(ncid,'elem_dim',varid)
      if(rcode.ne.0)return

      CALL NCVINQ(NCID,varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 230 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 230  CONTINUE
      RCODE= NF_GET_VARA_TEXT(NCID,varid,START,COUNT,elem_dim)
C
C    statements to fill line_dim                       
C
       rcode=NF_INQ_VARID(ncid,'line_dim',varid)
      if(rcode.ne.0)return

      CALL NCVINQ(NCID,varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 240 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 240  CONTINUE
      RCODE= NF_GET_VARA_TEXT(NCID,varid,START,COUNT,line_dim)
C
C    statements to fill imc                            
C
       rcode=NF_INQ_VARID(ncid,'imc',varid)
      if(rcode.ne.0)return

      CALL NCVINQ(NCID,varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 250 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 250  CONTINUE
      RCODE= NF_GET_VARA_TEXT(NCID,varid,START,COUNT,imc)
C
C    statements to fill eastWestCycles                 
C
       rcode=NF_INQ_VARID(ncid,'eastWestCycles',varid)
      if(rcode.ne.0)return

      CALL NCVINQ(NCID,varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 390 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 390  CONTINUE
      RCODE=NF_GET_VARA_INT(NCID,varid,START,COUNT,eastWestCycles)
C
C    statements to fill eastWestIncs                   
C
       rcode=NF_INQ_VARID(ncid,'eastWestIncs',varid)
      if(rcode.ne.0)return

      CALL NCVINQ(NCID,varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 400 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 400  CONTINUE
      RCODE=NF_GET_VARA_INT(NCID,varid,START,COUNT,eastWestIncs)
C
C    statements to fill northSouthCycles               
C
       rcode=NF_INQ_VARID(ncid,'northSouthCycles',varid)
      if(rcode.ne.0)return

      CALL NCVINQ(NCID,varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 410 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 410  CONTINUE
      RCODE=NF_GET_VARA_INT(NCID,varid,START,COUNT,northSouthCycles)
C
C    statements to fill northSouthIncs                 
C
       rcode=NF_INQ_VARID(ncid,'northSouthIncs',varid)
      if(rcode.ne.0)return

      CALL NCVINQ(NCID,varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 420 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 420  CONTINUE
      RCODE=NF_GET_VARA_INT(NCID,varid,START,COUNT,northSouthIncs)
C
C    statements to fill orbitAttitude                  
C
       rcode=NF_INQ_VARID(ncid,'orbitAttitude',varid)
      if(rcode.ne.0)return

      CALL NCVINQ(NCID,varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 430 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 430  CONTINUE
      RCODE=NF_GET_VARA_DOUBLE(NCID,varid,START,COUNT,orbitAttitude)
C
C    statements to fill Nx                             
C
       rcode=NF_INQ_VARID(ncid,'Nx',varid)
      if(rcode.ne.0)return

      CALL NCVINQ(NCID,varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 510 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 510  CONTINUE
      RCODE=NF_GET_VARA_INT(NCID,varid,START,COUNT,Nx)
      Nx4 = Nx(1)
C
C    statements to fill Ny                             
C
       rcode=NF_INQ_VARID(ncid,'Ny',varid)
      if(rcode.ne.0)return

      CALL NCVINQ(NCID,varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 520 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 520  CONTINUE
      RCODE=NF_GET_VARA_INT(NCID,varid,START,COUNT,Ny)
      Ny4=Ny(1)
C
C    statements to fill Lap                            
C
       rcode=NF_INQ_VARID(ncid,'Lap',varid)
      if(rcode.ne.0)return

      CALL NCVINQ(NCID,varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 530 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 530  CONTINUE
      RCODE=NF_GET_VARA_REAL(NCID,varid,START,COUNT,Lap)
C
C    statements to fill Lop                            
C
       rcode=NF_INQ_VARID(ncid,'Lop',varid)
      if(rcode.ne.0)return

      CALL NCVINQ(NCID,varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 540 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 540  CONTINUE
      RCODE=NF_GET_VARA_REAL(NCID,varid,START,COUNT,Lop)
C
C    statements to fill x_resolution                     
C
       rcode=NF_INQ_VARID(ncid,'x_resolution',varid)
      if(rcode.ne.0)return

      CALL NCVINQ(NCID,varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 550 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 550  CONTINUE
      RCODE=NF_GET_VARA_INT(NCID,varid,START,COUNT,x_resolution)
C
C    statements to fill y_resolution
C
       rcode=NF_INQ_VARID(ncid,'y_resolution',varid)
      if(rcode.ne.0)return

      CALL NCVINQ(NCID,varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 551 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 551  CONTINUE
      RCODE=NF_GET_VARA_INT(NCID,varid,START,COUNT,y_resolution)
C
C    statements to fill dx                             
C
       rcode=NF_INQ_VARID(ncid,'dx',varid)
      if(rcode.ne.0)return

      CALL NCVINQ(NCID,varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 580 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 580  CONTINUE
      RCODE=NF_GET_VARA_REAL(NCID,varid,START,COUNT,dx)
C
C    statements to fill dy                             
C
       rcode=NF_INQ_VARID(ncid,'dy',varid)
      if(rcode.ne.0)return

      CALL NCVINQ(NCID,varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 590 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 590  CONTINUE
      RCODE=NF_GET_VARA_REAL(NCID,varid,START,COUNT,dy)
C
C    statements to fill x_step
C
      rcode=NF_INQ_VARID(ncid,'x_step',varid)
      if(rcode.ne.0)return

      CALL NCVINQ(NCID,varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 680 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 680  CONTINUE
      RCODE=NF_GET_VARA_INT(NCID,varid,START,COUNT,x_step)
C
C    statements to fill y_step 
C
      rcode=NF_INQ_VARID(ncid,'y_step',varid)
      if(rcode.ne.0)return

      CALL NCVINQ(NCID,varid,DUMMY,NTP,NVDIM,VDIMS,NVS,RCODE)
      LENSTR=1
      DO 690 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 690  CONTINUE
      RCODE=NF_GET_VARA_INT(NCID,varid,START,COUNT,y_step)

      RETURN
      END
