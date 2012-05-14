      subroutine get_cdf_dims(nf_fid,record,x,y,istatus)

      include 'netcdf.inc'
      integer record,nf_fid
      integer istatus
      integer nf_status
      integer x, y

      istatus = 1
C
C  Fill all dimension values
C
C
C Get size of record
C
      nf_status = NF_INQ_DIMID(nf_fid,'record',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim record'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,record)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim record'
      endif
C
C Get size of x
C
      nf_status = NF_INQ_DIMID(nf_fid,'x',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim x'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,x)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim x'
      endif
C
C Get size of y
C
      nf_status = NF_INQ_DIMID(nf_fid,'y',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim y'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,y)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim y'
      endif

      istatus = 0

      return
      end

c ------------------------------------------------------------------

      subroutine read_sat_netcdf(nf_fid,record,x,y,image,istatus)

      include 'netcdf.inc'
      integer nf_fid,record,x,y
      integer image(x,y,record)
      integer istatus
c
c ------------------------------------------------------------
C
C     Variable        NETCDF Long Name
C      image        "GOES-8 satellite visible image"

      istatus = 1
 
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

      istatus = 0
 
      return
      end
c======================================================================
C
C  Subroutine to read the file "Image remapped to 5km FSL Conus C image " 
C
      subroutine read_netcdf_sat_head(nf_fid,record,
     +  Nx, Ny, center_id,process_id, wmo_sat_id,Dx,Dy,
     +  La1, Latin, Lo1, Lov, reftime, valtime, earth_shape,
     +  grid_name, grid_type, origin_name, process_name,
     +  wavelength, x_dim, y_dim, istatus)
C
      include 'netcdf.inc'
      integer record
      integer nf_fid, nf_vid
      integer nf_status,istatus
      integer Nx, Ny, center_id, process_id,
     +     wmo_sat_id(record)
      real Dx, Dy, La1, Latin, Lo1, Lov
      double precision reftime(record), valtime(record)
      character*132 origin_name
      character*132 x_dim
      character*132 y_dim
      character*132 earth_shape
      character*132 wavelength(record)
      character*132 grid_name
      character*132 process_name
      character*132 grid_type


      istatus = 1
C
C   Variables of type REAL
C
C     Variable        NETCDF Long Name
C      Dx           "x grid increment"
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
C
C     Variable        NETCDF Long Name
C      Dy           "y grid increment"
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
C
C     Variable        NETCDF Long Name
C      La1          "first latitude"
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
C      Latin        "conic tangent latitude"
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
C      Lo1          "first longitude"
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
C      Lov          "orientation of grid"
C
      nf_status = NF_INQ_VARID(nf_fid,'Lov',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var Lov'
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,Lov)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var Lov'
      endif

C   Variables of type INT
C
C
C     Variable        NETCDF Long Name
C      Nx           "number of x point"
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
C      Ny           "number of y points"
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
C
C     Variable        NETCDF Long Name
C      center_id    "center ID"
C
      nf_status = NF_INQ_VARID(nf_fid,'center_id',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var center_id'
      endif
      nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,center_id)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var center_id'
      endif
c
c ------------------------------------------------------------
C
C     Variable        NETCDF Long Name
C      image        "GOES-8 satellite visible image"
C
c     nf_status = NF_INQ_VARID(nf_fid,'image',nf_vid)
c    if(nf_status.ne.NF_NOERR) then
c       print *, NF_STRERROR(nf_status)
c       print *,'in var image'
c     endif
c     nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,image)
c     if(nf_status.ne.NF_NOERR) then
c       print *, NF_STRERROR(nf_status)
c       print *,'in var image'
c     endif
C --------------------------------------------------------------
C
C     Variable        NETCDF Long Name
C      process_id   "process ID"

      nf_status = NF_INQ_VARID(nf_fid,'process_id',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var process_id'
      endif
      nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,process_id)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var process_id'
      endif
C
C     Variable        NETCDF Long Name
C      wmo_sat_id   "WMO satellite id number"
C
      nf_status = NF_INQ_VARID(nf_fid,'wmo_sat_id',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var wmo_sat_id'
      else
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,wmo_sat_id)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var wmo_sat_id'
        endif
      endif

C   Variables of type DOUBLE
C
C
C     Variable        NETCDF Long Name
C      reftime      "reference time"
C
      nf_status = NF_INQ_VARID(nf_fid,'reftime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var reftime'
      endif
      nf_status = NF_GET_VAR_DOUBLE(nf_fid,nf_vid,reftime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var reftime'
      endif
C
C     Variable        NETCDF Long Name
C      valtime      "valid time"
C
      nf_status = NF_INQ_VARID(nf_fid,'valtime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var valtime'
      endif
      nf_status = NF_GET_VAR_DOUBLE(nf_fid,nf_vid,valtime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var valtime'
      endif


C   Variables of type CHAR
C
C
C     Variable        NETCDF Long Name
C      earth_shape  
C
      nf_status = NF_INQ_VARID(nf_fid,'earth_shape',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var earth_shape'
      endif
      nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,earth_shape)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var earth_shape'
      endif
C
C     Variable        NETCDF Long Name
C      grid_name    
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
C      grid_type    "GRIB-1 grid type"
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
C      origin_name  
C
      nf_status = NF_INQ_VARID(nf_fid,'origin_name',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var origin_name'
      endif
      nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,origin_name)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var origin_name'
      endif
C
C     Variable        NETCDF Long Name
C      process_name 
C
      nf_status = NF_INQ_VARID(nf_fid,'process_name',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var process_name'
      endif
      nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,process_name)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var process_name'
      endif
C
C     Variable        NETCDF Long Name
C      wavelength   "wavelength of satellite data"
C
      nf_status = NF_INQ_VARID(nf_fid,'wavelength',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var wavelength'
      endif
      nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,wavelength)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var wavelength'
      endif
C
C     Variable        NETCDF Long Name
C      x_dim        "longitude dimension"
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
C      y_dim        "latitude dimension"
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

      nf_status = nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
      endif

      return
      end
