cdis    Forecast Systems Laboratory
cdis    NOAA/OAR/ERL/FSL
cdis    325 Broadway
cdis    Boulder, CO     80303
cdis 
cdis    Forecast Research Division
cdis    Local Analysis and Prediction Branch
cdis    LAPS 
cdis 
cdis    This software and its documentation are in the public domain and 
cdis    are furnished "as is."  The United States government, its 
cdis    instrumentalities, officers, employees, and agents make no 
cdis    warranty, express or implied, as to the usefulness of the software 
cdis    and documentation for any purpose.  They assume no responsibility 
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis    
cdis    Permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  All modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making 
cdis    the modifications.  If significant modifications or enhancements 
cdis    are made to this software, the FSL Software Policy Manager  
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
      subroutine readcdf(csat_id,csat_type,chtype, Nx, Ny, 
     1n_elem,n_lines,r4_image, La1, Lo1,Dx, Dy, Latin,
     1LoV, validTime , data_file_name, istatus)

c     written by Dan Birkenheuer February 1994
c     J Smart   5/95    Included code to deal with gvar /public files
c     J Smart   6/96    Included code to deal with WFO satellite netcdf files
c     J Smart   3/97    Included subroutine rdblock_line_elem to read only the
c			appropriate sector of visible satellite data. Argument
c                       list lets subroutine know whether logical or I*2 data type.
c     J Smart   3/97    Converted _vis routine to be the only netCDF reader. Works
c			for vis, ir, wv and sounder. Block reading allows this.
c     J Smart   4/97    Adapted code further to work with gvar. NetCDF headers for
c                       raw gvar are different than for fsl-conus (ie., no La1, Lo1,
c                       Lov, or Latin).
c================================================================
      Include 'netcdf.inc'

      Integer n_elem,n_lines
      INTEGER RCODE
C
      real*4      r4_image(n_elem,n_lines)
      integer  ib
      integer       bi (2)
      equivalence (ib,bi (1) )
      Integer ndsize_x
      Integer ndsize_y,ndsize_ch
      character*200  data_file_name
      character*30 c_valtime
      character*30 c_Lov
      character*30 c_xres
      character*30 c_yres
      Character*3  csat_type
      character*3  chtype
      character*6  csat_id
      integer istatus
      INTEGER   validTime 
      Integer   ispec,istat
      Integer   iDx
      Integer   iDy
      INTEGER   Nx                             (   1)
      INTEGER   Ny                             (   1)
      REAL*4      La1                            (   1)
      REAL*4      Lo1                            (   1)
      REAL*4      Dx                             (   1)
      REAL*4      Dy                             (   1)
      REAL*4      Latin                          (   1)
      REAL*4      LoV                            (   1)
      INTEGER START(10)
      INTEGER COUNT(10)
      INTEGER VDIMS(10) !ALLOW UP TO 10 DIMENSIONS
      CHARACTER*31 DUMMY
      integer varid
c---------------------------------------------------------
c   code

      istatus = 0 ! bad istatus

      call NCPOPT(0)

      RCODE=NF_OPEN(data_file_name,NF_NOWRITE,NCID)
      if(rcode.ne.0) return

      rcode=NF_INQ_VARID(ncid,'image',varid)
      if(rcode.ne.0) return
C
C    statements to fill image                          
C
      Write(6,*)'Calling rdblock_line_elem - image read sub'

      Call rdblock_line_elem(csat_id,csat_type,chtype,
     &ncid,n_elem,n_lines,1,ndsize_x,ndsize_y,ndsize_ch,
     &r4_image,istatus)

      if(istatus .ne. 1)then
         write(6,*)'Error in rdblock_line_elem'
         return
      endif

      istatus = 0
c
c JSmart: 6-4-97.  WFO satellite netCDF files changed rather dramatically
c                  such that NO header info exists and we must now jump over
c                  the statements to read the header info. 
c
      if(csat_type.eq.'wfo')goto 201
C
C    statements to fill Nx
C
       rcode=NF_INQ_VARID(ncid,'Nx',varid)
      if(rcode.ne.0) return
      CALL NCVINQ(NCID, varid,'Nx',NTP,NVDIM,VDIMS,NVS,RCODE)
      if(rcode.ne.0) return

      LENSTR=1
      DO 150 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      if(rcode.ne.0) return

      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 150  CONTINUE
      RCODE=NF_GET_VARA_INT(NCID,varid,START,COUNT,Nx)
      if(rcode.ne.0) return

C    statements to fill Ny
C
       rcode=NF_INQ_VARID(ncid,'Ny',varid)
      if(rcode.ne.0) return
      CALL NCVINQ(NCID, varid,'Ny',NTP,NVDIM,VDIMS,NVS,RCODE)
      if(rcode.ne.0) return

      LENSTR=1
      DO 160 J=1,NVDIM
      CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
      if(rcode.ne.0) return

      LENSTR=LENSTR*NDSIZE
      START(J)=1
      COUNT(J)=NDSIZE
 160  CONTINUE
      RCODE=NF_GET_VARA_INT(NCID,varid,START,COUNT,Ny)
      if(rcode.ne.0) return
C
C    statements to fill validTime                      
C
c     CALL NCAINQ(ncid,varid,'valtime',itype,len,rcode)
c     if(rcode.ne.0)then
c        call NCAGTC(ncid,varid,'valtime',itype,ilen,rcode)
c        call NCVGT(ncid,varid,start,count,validTime,rcode)
c     else
c        call NCAINQ(ncid,varid,'validTime',itype,ilen,rcode)
c        call NCVGT(ncid,varid,start,count,validTime,rcode)
c     endif
c
c   statements to fill valtime
c
      c_valtime='valtime'                            !this for gvar/fsl-conus
      nvt = index(c_valtime,' ')-1
       rcode=NF_INQ_VARID(ncid,c_valtime(1:nvt),varid)
      if(rcode.ne.0)then
         if(csat_type.ne.'gvr')then
            c_valtime='validTime'                       !this for wfo
            nvt = index(c_valtime,' ')-1
            rcode=NF_INQ_VARID(ncid,c_valtime(1:nvt),varid)
            if(rcode .ne. 0)return
         else
            c_valtime='frameStartTime'                  !this for gvar/raw
            nvt = index(c_valtime,' ')-1
            rcode=NF_INQ_VARID(ncid,c_valtime(1:nvt),varid)
            if(rcode .ne. 0)return
         endif
      endif

      CALL NCVINQ(NCID, varid,c_valtime,NTP,NVDIM,VDIMS,NVS,RCODE)

      if(rcode.ne.0) return

      LENSTR=1
      RCODE=NF_GET_VAR_INT(NCID,varid,validtime)
      if(rcode.ne.0) return
C
C    statements to fill Dx
C
      if(csat_type.eq.'cdf'.or.csat_type.eq.'wfo')then
         c_xres = 'Dx'
         c_yres = 'Dy'
      elseif(csat_type.eq.'gvr')then
         c_xres = 'x_resolution'
         c_yres = 'y_resolution'
      endif
      nxr=index(c_xres,' ')
      nyr=index(c_yres,' ')

       rcode=NF_INQ_VARID(ncid,c_xres(1:nxr),varid)
      if(rcode.ne.0) return
      CALL NCVINQ(NCID, varid,c_xres(1:nxr),NTP,NVDIM,VDIMS,NVS,RCODE)

      if(rcode.ne.0) return

      LENSTR=1
      DO 210 J=1,NVDIM
         CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
         if(rcode.ne.0) return

         LENSTR=LENSTR*NDSIZE
         START(J)=1
         COUNT(J)=NDSIZE
210   CONTINUE
      if(csat_type.eq.'cdf'.or.csat_type.eq.'wfo')then

      RCODE=NF_GET_VARA_REAL(NCID,varid,START,COUNT,Dx)
         if(rcode.ne.0) return

      elseif(csat_type.eq.'gvr')then

      RCODE=NF_GET_VARA_INT(NCID,varid,START,COUNT,iDx)
         if(rcode.ne.0) return

       Dx(1) = float(iDx)

      endif

C
C    statements to fill Dy
C
       rcode=NF_INQ_VARID(ncid,c_yres(1:nyr),varid)
      if(rcode.ne.0) return
      CALL NCVINQ(NCID, varid,c_yres(1:nyr),NTP,NVDIM,VDIMS,NVS,RCODE)
      if(rcode.ne.0) return

      LENSTR=1
      DO 220 J=1,NVDIM
         CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
         if(rcode.ne.0) return

         LENSTR=LENSTR*NDSIZE
         START(J)=1
         COUNT(J)=NDSIZE
220   CONTINUE
      if(csat_type.eq.'cdf'.or.csat_type.eq.'wfo')then

      RCODE=NF_GET_VARA_REAL(NCID,varid,START,COUNT,Dy)
         if(rcode.ne.0) return

      elseif(csat_type .eq. 'gvr')then

      RCODE=NF_GET_VARA_INT(NCID,varid,START,COUNT,iDy)
         if(rcode.ne.0) return

       Dy(1) = float(iDy)

      endif
C
C    statements to fill Latin
C
      if(csat_type.eq.'cdf'.or.csat_type.eq.'wfo')then
          rcode=NF_INQ_VARID(ncid,'Latin',varid)
         if(rcode.ne.0) return
         CALL NCVINQ(NCID, varid,'Latin',NTP,NVDIM,VDIMS,NVS,RCODE)
         if(rcode.ne.0) return
         LENSTR=1
         DO 230 J=1,NVDIM
            CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
            if(rcode.ne.0) return

            LENSTR=LENSTR*NDSIZE
            START(J)=1
            COUNT(J)=NDSIZE
230      CONTINUE
         RCODE=NF_GET_VARA_REAL(NCID,varid,START,COUNT,Latin)
         if(rcode.ne.0) return
C
C    statements to fill LoV
C
          rcode=NF_INQ_VARID(ncid,'LoV',varid)
         c_Lov = 'LoV'
         if(rcode.ne.0)then
            rcode=NF_INQ_VARID(ncid,'Lov',varid)
            c_Lov = 'Lov'
            if(rcode.ne.0)then
               write(6,*)'Error getting variable LoV/Lov'
               goto 241
            endif
         endif
 
         CALL NCVINQ(NCID, varid,c_Lov,NTP,NVDIM,VDIMS,NVS,RCODE)
         if(rcode.ne.0) return

         LENSTR=1
         DO 240 J=1,NVDIM
            CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
            if(rcode.ne.0) return

            LENSTR=LENSTR*NDSIZE
            START(J)=1
            COUNT(J)=NDSIZE
240      CONTINUE
            RCODE=NF_GET_VARA_REAL(NCID,varid,START,COUNT,LoV)
         if(rcode.ne.0) return
C
C    statements to fill La1
C
241      continue

          rcode=NF_INQ_VARID(ncid,'La1',varid)
         if(rcode.ne.0) return
         CALL NCVINQ(NCID, varid,'La1',NTP,NVDIM,VDIMS,NVS,RCODE)
         if(rcode.ne.0) return

         LENSTR=1
         DO 170 J=1,NVDIM
            CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
            if(rcode.ne.0) return

            LENSTR=LENSTR*NDSIZE
            START(J)=1
            COUNT(J)=NDSIZE
170      CONTINUE
         RCODE=NF_GET_VARA_REAL(NCID,varid,START,COUNT,La1)
         if(rcode.ne.0) return
C
C    statements to fill Lo1
C
          rcode=NF_INQ_VARID(ncid,'Lo1',varid)
         if(rcode.ne.0) return
         CALL NCVINQ(NCID, varid,'Lo1',NTP,NVDIM,VDIMS,NVS,RCODE)
         if(rcode.ne.0) return

         LENSTR=1
         DO 180 J=1,NVDIM
            CALL NCDINQ(NCID,VDIMS(J),DUMMY,NDSIZE,RCODE)
            if(rcode.ne.0) return

            LENSTR=LENSTR*NDSIZE
            START(J)=1
            COUNT(J)=NDSIZE
180      CONTINUE
         RCODE=NF_GET_VARA_REAL(NCID,varid,START,COUNT,Lo1)
         if(rcode.ne.0) return
C
      endif !test for gvr

201   istatus = 1  ! ok!

      rcode= NF_CLOSE(ncid)

      Return
      END
