      subroutine rd_static_attr_sub(staticfile, Nx, Ny
     .,La1, Latin1, Latin2, Lo1, LoV, Dx, Dy
     .,c8_maproj,istatus)
C
      include 'netcdf.inc'
      integer nf_fid, nf_vid, nf_status
      integer Nx, Ny
      integer istatus
      real Dx, Dy, La1, Latin1, Latin2, Lo1, LoV

      character*132 grid_type
      character*132 grid_name
      character*132 earth_shape
      character*8   c8_maproj
      character*(*) staticfile
C
C  Open netcdf File for reading
C
      istatus = 0
 
      nf_status = NF_OPEN(TRIM(staticfile),NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'NF_OPEN ',TRIM(staticfile)
        return
      endif
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
        return
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,Dx)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var Dx'
        return
      endif
C
C     Variable        NETCDF Long Name
C      Dy           "y grid increment"
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
        print *,'in var Dy'
        return
      endif
C
C     Variable        NETCDF Long Name
C      La1          "first latitude"
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
        print *,'in var La1'
        return
      endif
C
C     Variable        NETCDF Long Name
C      Latin1       "orientation of grid"
C
      nf_status = NF_INQ_VARID(nf_fid,'Latin1',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var Latin1'
        return
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,Latin1)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var Latin1'
        return
      endif
C
C     Variable        NETCDF Long Name
C      Latin2       "orientation of grid"
C
      nf_status = NF_INQ_VARID(nf_fid,'Latin2',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var Latin2'
        return
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,Latin2)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var Latin2'
        return
      endif
C
C     Variable        NETCDF Long Name
C      Lo1          "first longitude"
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
        print *,'in var Lo1'
        return
      endif
C
C     Variable        NETCDF Long Name
C      LoV          "orientation of grid"
C
      nf_status = NF_INQ_VARID(nf_fid,'LoV',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var LoV'
        return
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,LoV)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var LoV'
        return
      endif
C
C     Variable        NETCDF Long Name
C      grid_spacing 
C
      nf_status = NF_INQ_VARID(nf_fid,'grid_spacing',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var grid_spacing'
        return
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,grid_spacing)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var grid_spacing'
        return
      endif
C
C   Variables of type INT
C
C
C     Variable        NETCDF Long Name
C      Nx           "number of x points"
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
        print *,'in var Nx'
        return
      endif
C
C     Variable        NETCDF Long Name
C      Ny           "number of y points"
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
        print *,'in var Ny'
        return
      endif
C
C      earth_shape  
C
      nf_status = NF_INQ_VARID(nf_fid,'earth_shape',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var earth_shape'
        return
      endif
      nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,earth_shape)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var earth_shape'
        return 
      endif
C
C     Variable        NETCDF Long Name
C      grid_name    
C
      nf_status = NF_INQ_VARID(nf_fid,'grid_name',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var grid_name'
        return
      endif
      nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,grid_name)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var grid_name'
        return
      endif
C
C     Variable        NETCDF Long Name
C      grid_type    "GRIB-1 grid type"
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
        print *,'in var grid_type'
        return
      endif

      call downcase(grid_type,grid_type)
      call s_len(grid_type,lgt)
      if(grid_type(1:lgt).eq.'tangential'.or.
     &grid_type(1:lgt).eq.'secant')then
         c8_maproj='lambert'
      elseif(grid_type(1:lgt).eq.'mercator')then
         c8_maproj='mercator'
      else
         c8_maproj='polar'
      endif
      call s_len(grid_type,lgt)

      if(TRIM(grid_type).eq.'tangential lambert conformal'
     &.or.TRIM(grid_type).eq.'secant lambert conformal')then
         c8_maproj='lambert'
      elseif(TRIM(grid_type).eq.'mercator')then
         c8_maproj='mercator'
      elseif(TRIM(grid_type).eq.'polar')then
         c8_maproj='polar'
      endif
C
      nf_status = nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
        return
      endif

      istatus = 1
      return
      end
