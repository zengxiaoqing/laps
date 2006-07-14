      subroutine read_nowrad_cdf(ctype,cfname_in,lines,elems,
     + dlat,dlon,La1,Lo1,La2,Lo2,centerLon,topLat,validtime,
     + Dx,Dy,Lov, Latin, image, istatus)

      include 'netcdf.inc'
      integer elems, lines,nf_fid, nf_vid, nf_status
      integer Nx, Ny, image( elems, lines), validTime
      real Dx, Dy, La1, La2, Lo1, Lo2, centerLon, diffLon,
     +     dlat, dlon, topLat, LoV, Latin

      character*200 x_dim
      character*200 y_dim
      character*200 grid_name
      character*200 grid_type
      character*(*)  cfname_in
      character*(*)  ctype
C
      call read_nowrad_head(ctype,cfname_in,Nx,Ny,validTime,
     + Dx, Dy, La1, La2, Lo1, Lo2, centerLon, diffLon,dlat,dlon,
     + topLat, grid_name, grid_type,x_dim,y_dim,LoV,Latin)

      call read_wsi_image(cfname_in,elems,lines,image)

      return
      end
C
C  ===============================================================
C  Subroutine to read the file "WSI NOWRAD Data Header" 
C  ===============================================================
C
      subroutine read_nowrad_head(ctype,cfname_in,Nx,Ny, 
     +     validTime, Dx, Dy, La1, La2, Lo1, Lo2, centerLon, diffLon, 
     +     radsPerElem, radsPerLine, topLat, grid_name, grid_type, 
     +     x_dim, y_dim, LoV, Latin)
C
      include 'netcdf.inc'
      integer elems, lines,nf_fid, nf_vid, nf_status
      integer Nx, Ny, validTime
      real Dx, Dy, La1, La2, Lo1, Lo2, centerLon, diffLon,
     +     radsPerElem, radsPerLine, topLat, LoV, Latin

      real*8 valTime
      integer attlen
      character c_atvalue*80
      character*(*) x_dim
      character*(*) y_dim
      character*(*) grid_name
      character*(*) grid_type
      character*(*) cfname_in
      character*(*) ctype
C
C  Open netcdf File for reading
C
      call s_len(cfname_in,nc)
      nf_status = NF_OPEN(cfname_in,NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'NF_OPEN ', cfname_in(1:nc)
        return
      endif

C   Variables of type REAL
C
C     Variable        NETCDF Long Name
C      Dx           
C
        nf_status = NF_INQ_VARID(nf_fid,'Dx',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var Dx'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,Dx)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var Dx'
      endif

      call NCAINQ(nf_fid,nf_vid,'units',itype,attlen,istatus)
      call NCAGTC(nf_fid,nf_vid,'units',c_atvalue,attlen,istatus)
      if(istatus.ne.0)then
         write(6,*)'Error getting attribute - Dx'
      endif

C
C     Variable        NETCDF Long Name
C      Dy           
C
        nf_status = NF_INQ_VARID(nf_fid,'Dy',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var Dy'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,Dy)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var Dy'
      endif
c
c wsi-type is degrees. wfo-type may be kilometers. Dx and Dy must = m.
c
      if(c_atvalue(1:4).eq.'kilo')then
         Dx=Dx*1000.
         Dy=Dy*1000.
      elseif(c_atvalue(1:4).eq.'degr')then
         Dx=Dx*111.1*1000.
         Dy=Dy*111.1*1000.
      endif
C
C     Variable        NETCDF Long Name
C      La1          
C
        nf_status = NF_INQ_VARID(nf_fid,'La1',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var La1'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,La1)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var La1'
      endif
C
C     Variable        NETCDF Long Name
C      La2          
C
        nf_status = NF_INQ_VARID(nf_fid,'La2',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var La2'
      else
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,La2)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var La2'
        endif
      endif
C
C     Variable        NETCDF Long Name
C      Lo1          
C
        nf_status = NF_INQ_VARID(nf_fid,'Lo1',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var Lo1'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,Lo1)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var Lo1'
      endif
C
C     Variable        NETCDF Long Name
C      Lo2          
C
        nf_status = NF_INQ_VARID(nf_fid,'Lo2',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var Lo2'
      else
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,Lo2)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var Lo2'
        endif
      endif

c ******** wsi type variables ***********

      if(ctype.eq.'wsi')then

C
C     Variable        NETCDF Long Name
C      centerLon    "Center Longitude"
C
          nf_status = NF_INQ_VARID(nf_fid,'centerLon',nf_vid)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var centerLon'
        endif
          nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,centerLon)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var centerLon'
        endif
C
C     Variable        NETCDF Long Name
C      diffLon      "Difference Longitude"
C
          nf_status = NF_INQ_VARID(nf_fid,'diffLon',nf_vid)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var diffLon'
        endif
          nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,diffLon)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var diffLon'
        endif
C
C     Variable        NETCDF Long Name
C      radsPerElem  "Radians per element"
C
          nf_status = NF_INQ_VARID(nf_fid,'radsPerElem',nf_vid)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var radsPerElem'
        endif
          nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,radsPerElem)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var radsPerElem'
        endif
C
C     Variable        NETCDF Long Name
C      radsPerLine  "Radians per line"
C
          nf_status = NF_INQ_VARID(nf_fid,'radsPerLine',nf_vid)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var radsPerLine'
        endif
          nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,radsPerLine)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var radsPerLine'
        endif
C
C     Variable        NETCDF Long Name
C      topLat       "Top Latitude"
C
          nf_status = NF_INQ_VARID(nf_fid,'topLat',nf_vid)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var topLat'
        endif
          nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,topLat)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var topLat'
        endif

C ------------------------------

      elseif(ctype.eq.'wfo')then

C
C     Variable        NETCDF Long Name
C      Latin       "Latitude of Tangent"
C
          nf_status = NF_INQ_VARID(nf_fid,'Latin',nf_vid)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var Latin'
        endif
          nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,Latin)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var Latin'
        endif
C
C     Variable        NETCDF Long Name
C      LoV       "Longitude of Vertical"
C
          nf_status = NF_INQ_VARID(nf_fid,'LoV',nf_vid)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var LoV'
        endif
          nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,LoV)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var LoV'
        endif

      else

        print*,'Invalid radar type in ctype ',ctype

      endif

C   Variables of type INT
C
C
C     Variable        NETCDF Long Name
C      Nx           "Number of x Points"
C
        nf_status = NF_INQ_VARID(nf_fid,'Nx',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var Nx'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,Nx)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var Nx'
      endif
C
C     Variable        NETCDF Long Name
C      Ny           "Number of y Points"
C
        nf_status = NF_INQ_VARID(nf_fid,'Ny',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var Ny'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,Ny)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var Ny'
      endif

      if(ctype.eq.'wsi')then
C
C     Variable        NETCDF Long Name
C      validTime    "Valid Time"
C
          nf_status = NF_INQ_VARID(nf_fid,'validTime',nf_vid)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var validTime'
        endif
          nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,validTime)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var validTime'
        endif

      elseif(ctype.eq.'wfo')then

C
C     Variable        NETCDF Long Name
C      validTime    "Valid Time"
C
        nf_status=NF_INQ_VARID(nf_fid,'valtime',nf_vid)
        if(nf_status.ne.0)then
           nf_status=NF_INQ_VARID(nf_fid,'validTime',nf_vid)
           if(nf_status.ne.0)then
              write(6,*)'Error getting variable - valtime'
              return
           endif
        endif
        nf_status=NF_GET_VAR1_DOUBLE(nf_fid,nf_vid,1,valtime)
        validTime=int(valTime)

      endif

C   Variables of type CHAR
C
C
C     Variable        NETCDF Long Name
C      grid_name    "Grid Name"
C
        nf_status = NF_INQ_VARID(nf_fid,'grid_name',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var grid_name'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,grid_name)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var grid_name'
      endif
C
C     Variable        NETCDF Long Name
C      grid_type    "Grid Type"
C
        nf_status = NF_INQ_VARID(nf_fid,'grid_type',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var grid_type'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,grid_type)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var grid_type'
      endif
C
C     Variable        NETCDF Long Name
C      x_dim        "Longitude Dimension"
C
        nf_status = NF_INQ_VARID(nf_fid,'x_dim',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var x_dim'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,x_dim)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var x_dim'
      endif
C
C     Variable        NETCDF Long Name
C      y_dim        "Latitude Dimension"
C
        nf_status = NF_INQ_VARID(nf_fid,'y_dim',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var y_dim'
      endif
        nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,y_dim)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var y_dim'
      endif
c
c close file
c
      nf_status = nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
      endif

      return
      end
C
C ================================================================
C
      subroutine read_wsi_image(cfname_in,elems,lines,image)

      include 'netcdf.inc'

      integer elems, lines,nf_fid, nf_vid, nf_status
      integer image(elems,lines)
      character*(*) cfname_in
C
C  Open netcdf File for reading
C
      call s_len(cfname_in,nc)
      nf_status = NF_OPEN(cfname_in,NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'NF_OPEN ', cfname_in(1:nc)
        return
      endif
C
C     Variable        NETCDF Long Name
C      image        "Image Pixel Values"
C
        nf_status = NF_INQ_VARID(nf_fid,'image',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var image'
      endif
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,image)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var image'
      endif

c close file
c
      nf_status = nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
      endif

      return
      end
c
c ============================================================
c Routine to convert NETCDF data base (bytes) to NOWRAD scaled
c integers.
c ============================================================

      subroutine cvt_wsi_nowrad(ctype,elems,lines,image
     +,nsites_present,present_site_loc_i,present_site_loc_j
     +,nsites_absent,absent_site_loc_i,absent_site_loc_j
     +,max_sites,istatus)
c
c
      implicit none

      integer elems,lines
      integer image(elems,lines)

      integer max_sites

      integer nsites_present
      real    present_site_loc_i(max_sites)
      real    present_site_loc_j(max_sites)

      integer nsites_absent
      real    absent_site_loc_i(max_sites)
      real    absent_site_loc_j(max_sites)

      integer istatus
      integer imax_image_value,imin_image_value,i,j
      integer icount_out, icount_bad
      integer bad_data_flag
      parameter(bad_data_flag=255) 
 
      character ctype*(*)

      if(ctype.eq.'wfo')then
         do j=1,lines
         do i=1,elems
            if(image(i,j).lt.0)image(i,j)=256+image(i,j)
            image(i,j)=image(i,j)/16
         enddo
         enddo
      endif
 
      imax_image_value = 0
      imin_image_value = 255
      icount_bad=0
      icount_out=0
      do j=1,lines
         do i=1,elems
            if(image(i,j) .gt. imax_image_value)
     +           imax_image_value = image(i,j)
            if(image(i,j) .lt. imin_image_value)
     +           imin_image_value = image(i,j)
            if(image(i,j) .ne. bad_data_flag)then
               if((image(i,j).gt.64).or.(image(i,j).lt.0))then
                  icount_out = icount_out +1
                  image(i,j) = bad_data_flag
c     write(6,*) i, j, i_value
               endif

            else
               icount_bad=icount_bad+1
               image(i,j) = bad_data_flag
            endif
         enddo
      enddo
c
      write(6,*)'Number of bad data points (> ',bad_data_flag,' )'
      write(6,*)'prior to calling c_scan_adjust: ',icount_bad
      write(6,*)'Max value found in image array: ',imax_image_value
      write(6,*)'Min value found in image array: ',imin_image_value
      write(6,*)
      write(6,*)'Data found out-of-bounds (icount_out) ',icount_out
      if(icount_out.gt.0)then
         istatus = -1
         return
      endif


C       new C version of scan_adjust implemented 05-apr-96
ccc      call c_scan_adjust(image,lines,elems,bad_data_flag)

      nsites_present=0
      nsites_absent =0
      if(ctype.eq.'wsi')then
         do j=1,lines
         do i=1,elems
            if(image(i,j).ne.bad_data_flag)then
               if(image(i,j).gt.15.and.image(i,j).lt.32)then
                  nsites_present=nsites_present+1
                  present_site_loc_i(nsites_present)=float(i)
                  present_site_loc_j(nsites_present)=float(j)
               elseif(image(i,j).ge.32 .and.image(i,j).lt.48)then
                  nsites_absent=nsites_absent+1
                  absent_site_loc_i(nsites_absent)=float(i)
                  absent_site_loc_j(nsites_absent)=float(j)
               endif
               image(i,j)=mod(image(i,j),16)
            endif
         enddo
         enddo
      endif

 9999 return
      end
c
c ----------------------------------------------------------
c
      subroutine read_nowrad_dims(ctype,cfname,elems,lines)

      include 'netcdf.inc'

      integer elems,lines
      character*(*) cfname
      character*(*) ctype
      character*10  cvarname(2)

C
C  Open netcdf File for reading
C
      call s_len(cfname,nc)
      nf_status = NF_OPEN(cfname,NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'NF_OPEN ', cfname(1:nc)
        return
      endif
      cvarname(1)='x'
      cvarname(2)='y'
      if(ctype.eq.'wsi')then
         cvarname(1)='elems'
         cvarname(2)='lines'
      endif
      call s_len(cvarname(1),il)
C
C Get size of elems
C
      nf_status = NF_INQ_DIMID(nf_fid,cvarname(1)(1:il),nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim elems'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,elems)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim elems'
      endif
C
C Get size of lines
C
      nf_status = NF_INQ_DIMID(nf_fid,cvarname(2)(1:il),nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim lines'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,lines)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim lines'
      endif
c
c close file
c
      nf_status = nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
      endif

      return
      end

