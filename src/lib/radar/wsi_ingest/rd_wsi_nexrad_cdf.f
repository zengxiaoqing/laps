      subroutine rd_wsi_3dradar_cdf(cfname,lines,elems,
     &radsPerLine,radsPerElem,la1,lo1,la2,lo2,centerLon,
     &validTime,dataLevel,data_levels,num_levels,level_prefix,
     &image,istatus)

      include 'netcdf.inc'
c     integer dataLevel, elems, lines, nf_fid, nf_vid, nf_status
      integer dataLevel, elems, lines
      integer image(elems,lines)
      integer istatus
      integer validtime

      character*(*) cfname
      character  level_prefix(dataLevel)*50

      real Dx,Dy,La1,Lo1,centerLon
      real radsPerElem,radsPerLine
c    +   , topLat

      istatus = 1

      call rd_wsi3d_sub(cfname,dataLevel,elems,lines,
     .   dx,dy,la1,lo1,la2,lo2,centerlon,
     .   radsPerLine,radsPerElem,num_levels,data_levels,
     .   level_prefix,image,validTime,istatus)
      if(istatus.ne.0)then
         print*,'Error reading wsi 3d radar file'
         return
      endif

      istatus = 0

      return
      end
C
      subroutine rd_wsi3d_sub(cfname,dataLevel,elems,lines,
     .   dx,dy,la1,lo1,la2,lo2,centerLon,radsPerLine,
     .   radsPerElem,num_levels,data_levels,level_prefix,image,
     .   validTime,istatus)


      include 'netcdf.inc'
      integer dataLevel, elems, lines
      character*50 grid_name
      character*50 grid_type
      character*50 level_prefix(dataLevel)
      character*20 product_units
      character*50 x_dim
      character*50 y_dim
      character*(*) cfname
      integer Nx, Ny, data_levels(dataLevel),image(elems,lines), 
     +   num_levels, validTime
      real Dx,Dy,La1,La2,Lo1,Lo2,centerLon,diffLon,
     +   radsPerElem,  radsPerLine, topLat

      call read_cdf_wsi3d_head(cfname,dataLevel, Dx, Dy,
     +    La1, La2, Lo1, Lo2, Nx, Ny, centerLon,
     +    data_levels,
     +    diffLon, grid_name, grid_type, level_prefix,
     +    num_levels, product_units, radsPerElem, radsPerLine,
     +    topLat, validTime, x_dim, y_dim, istatus)
C

      call read_cdf_wsi3d_image(cfname,elems,lines,image,istatus)

      return
      end
C
C  Subroutine to read the file header "WSI NEXRAD Data" 
C
      subroutine read_cdf_wsi3d_head(cfname,dataLevel, Dx, Dy,
     +    La1, La2, Lo1, Lo2, Nx, Ny,  centerLon,
     +    data_levels,diffLon,grid_name,grid_type,level_prefix,
     +    num_levels, product_units, radsPerElem, radsPerLine,
     +    topLat, validTime, x_dim, y_dim,istatus)
      include 'netcdf.inc'
c     integer dataLevel, elems, lines, nf_fid, nf_vid, nf_status
      integer dataLevel, nf_fid, nf_vid, nf_status

      character*50 grid_name
      character*50 grid_type
      character*50 level_prefix(dataLevel)
      character*20 product_units
      character*50 x_dim
      character*50 y_dim
      character*(*) cfname
      integer Nx, Ny, data_levels(dataLevel), 
     +   num_levels, validTime
      real Dx, Dy, La1, La2, Lo1, Lo2, centerLon,
     +   diffLon, radsPerElem,  radsPerLine, topLat
C
C  Open netcdf File for reading
C
      nf_status = NF_OPEN(cfname,NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        call s_len(cfname,nf)
        print *, NF_STRERROR(nf_status)
        print *,'NF_OPEN ',cfname(1:nf)
        return
      endif
C
C     Variable        NETCDF Long Name
C      grid_name    "Grid Name" 
C
      istatus=1
      nf_status = NF_INQ_VARID(nf_fid,'grid_name',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var grid_name'
        return
      endif
      nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,grid_name)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ grid_name '
        return
      endif
C
C     Variable        NETCDF Long Name
C      grid_type    "Grid Type" 
C
      nf_status = NF_INQ_VARID(nf_fid,'grid_type',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var grid_type'
        return
      endif
      nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,grid_type)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ grid_type '
        return
      endif
C
C     Variable        NETCDF Long Name
C      level_prefix "Level Prefix" 
C
      nf_status = NF_INQ_VARID(nf_fid,'level_prefix',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var level_prefix'
        return
      endif
      nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,level_prefix)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ level_prefix '
        return
      endif
C
C     Variable        NETCDF Long Name
C      product_units"Product Units" 
C
      nf_status = NF_INQ_VARID(nf_fid,'product_units',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var product_units'
        return
      endif
      nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,product_units)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ product_units '
        return
      endif
C
C     Variable        NETCDF Long Name
C      x_dim        "Longitude Dimension" 
C
      nf_status = NF_INQ_VARID(nf_fid,'x_dim',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var x_dim'
        return
      endif
      nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,x_dim)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ x_dim '
        return
      endif
C
C     Variable        NETCDF Long Name
C      y_dim        "Latitude Dimension" 
C
      nf_status = NF_INQ_VARID(nf_fid,'y_dim',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var y_dim'
        return
      endif
      nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,y_dim)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ y_dim '
        return
      endif
C
C     Variable        NETCDF Long Name
C      Nx           "Number of x Points" 
C
      nf_status = NF_INQ_VARID(nf_fid,'Nx',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var Nx'
        return
      endif
      nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,Nx)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ Nx '
        return
      endif
C
C     Variable        NETCDF Long Name
C      Ny           "Number of y Points" 
C
      nf_status = NF_INQ_VARID(nf_fid,'Ny',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var Ny'
        return
      endif
      nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,Ny)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ Ny '
        return
      endif
C
C     Variable        NETCDF Long Name
C      data_levels  "Data Levels" 
C
      nf_status = NF_INQ_VARID(nf_fid,'data_levels',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var data_levels'
        return
      endif
      nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,data_levels)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ data_levels '
        return
      endif
C
C     Variable        NETCDF Long Name
C      num_levels   "Number of Levels" 
C
      nf_status = NF_INQ_VARID(nf_fid,'num_levels',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var num_levels'
        return
      endif
      nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,num_levels)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ num_levels '
        return
      endif
C
C     Variable        NETCDF Long Name
C      validTime    "Valid Time" 
C
      nf_status = NF_INQ_VARID(nf_fid,'validTime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var validTime'
        return
      endif
      nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,validTime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ validTime '
        return
      endif
C
C     Variable        NETCDF Long Name
C      Dx           
C
      nf_status = NF_INQ_VARID(nf_fid,'Dx',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var Dx'
        return
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,Dx)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ Dx '
        return
      endif
C
C     Variable        NETCDF Long Name
C      Dy           
C
      nf_status = NF_INQ_VARID(nf_fid,'Dy',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var Dy'
        return
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,Dy)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ Dy '
        return
      endif
C
C     Variable        NETCDF Long Name
C      La1          
C
      nf_status = NF_INQ_VARID(nf_fid,'La1',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var La1'
        return
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,La1)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ La1 '
        return
      endif
C
C     Variable        NETCDF Long Name
C      La2          
C
      nf_status = NF_INQ_VARID(nf_fid,'La2',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var La2'
        return
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,La2)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ La2 '
        return
      endif
C
C     Variable        NETCDF Long Name
C      Lo1          
C
      nf_status = NF_INQ_VARID(nf_fid,'Lo1',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var Lo1'
        return
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,Lo1)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ Lo1 '
        return
      endif
C
C     Variable        NETCDF Long Name
C      Lo2          
C
      nf_status = NF_INQ_VARID(nf_fid,'Lo2',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var Lo2'
        return
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,Lo2)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ Lo2 '
        return
      endif
C
C     Variable        NETCDF Long Name
C      centerLon    "Center Longitude" 
C
      nf_status = NF_INQ_VARID(nf_fid,'centerLon',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var centerLon'
        return
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,centerLon)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ centerLon '
        return
      endif
C
C     Variable        NETCDF Long Name
C      diffLon      "Difference Longitude" 
C
      nf_status = NF_INQ_VARID(nf_fid,'diffLon',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var diffLon'
        return
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,diffLon)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ diffLon '
        return
      endif
C
C     Variable        NETCDF Long Name
C      radsPerElem  "Radians per element" 
C
      nf_status = NF_INQ_VARID(nf_fid,'radsPerElem',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var radsPerElem'
        return
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,radsPerElem)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ radsPerElem '
        return
      endif
C
C     Variable        NETCDF Long Name
C      radsPerLine  "Radians per line" 
C
      nf_status = NF_INQ_VARID(nf_fid,'radsPerLine',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var radsPerLine'
        return
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,radsPerLine)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ radsPerLine '
        return
      endif
C
C     Variable        NETCDF Long Name
C      topLat       "Top Latitude" 
C
      nf_status = NF_INQ_VARID(nf_fid,'topLat',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var topLat'
        return
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,topLat)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in NF_GET_VAR_ topLat '
        return
      endif
      nf_status = nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
        return
      endif

      istatus = 0

      return
      end

c
c --------------------------------------------------------
c
      subroutine get_wsi_3drad_dims(cfname,dataLevel,elems,lines,
     &                              istatus)

      include 'netcdf.inc'
      integer dataLevel, elems, lines, nf_fid, nf_vid, nf_status
      integer istatus
      character*(*) cfname

      istatus=1  !return error status 
C
C  Open netcdf File for reading
C
      nf_status = NF_OPEN(cfname,NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        call s_len(cfname,nf)
        print *, NF_STRERROR(nf_status)
        print *,'NF_OPEN',cfname(1:nf)
        return
      endif
C
C  Fill all dimension values
C
C
C Get size of dataLevel
C
      nf_status = NF_INQ_DIMID(nf_fid,'dataLevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim dataLevel'
        return
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,dataLevel)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim dataLevel'
        return
      endif
C
C Get size of elems
C
      nf_status = NF_INQ_DIMID(nf_fid,'elems',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim elems'
        return
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,elems)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim elems'
        return
      endif
C
C Get size of lines
C
      nf_status = NF_INQ_DIMID(nf_fid,'lines',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim lines'
        return
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,lines)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim lines'
        return
      endif
c
c close file
c
      nf_status = nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
        return
      endif

      istatus = 0
      return
      end
c
c ----------------------------------------------------
c
      subroutine read_cdf_wsi3d_image(cfname,
     &              elems,lines,image,istatus)
C
C     Variable        NETCDF Long Name
C      image        "Image Pixel Values"
C
      include 'netcdf.inc'
      integer elems,lines
      integer image(elems,lines)
      character*(*) cfname
C
C  Open netcdf File for reading
C
      nf_status = NF_OPEN(cfname,NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        call s_len(cfname,nf)
        print *, NF_STRERROR(nf_status)
        print *,'NF_OPEN ',cfname(1:nf)
        return
      endif

      nf_status = NF_INQ_VARID(nf_fid,'image',nf_vid)
      if(nf_status.ne.NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *,'in var image'
        return
      endif
      nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,image)
      if(nf_status.ne.NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *,'in NF_GET_VAR_ image '
        return
      endif
c
c close file
c
      nf_status = nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
        return
      endif

      return
      end
