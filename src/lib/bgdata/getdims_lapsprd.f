      subroutine getdims_lapsprd(fullname,x,y,z,istatus)
      implicit none
      include 'netcdf.inc'
      integer record, x, y, z
      integer nf_fid, nf_vid, nf_status, ln
      integer istatus
      character*(*) fullname
C
C  Open netcdf File for reading
C
      istatus = 0

      call s_len(fullname,ln)
      nf_status = NF_OPEN(fullname,NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'NF_OPEN ',fullname(1:ln)
        return
      endif
C
C Get size of x
C
      nf_status = NF_INQ_DIMID(nf_fid,'x',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim x'
        return
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,x)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim x'
        return
      endif
C
C Get size of y
C
      nf_status = NF_INQ_DIMID(nf_fid,'y',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim y'
        return
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,y)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim y'
        return
      endif
C
C Get size of z
C
      nf_status = NF_INQ_DIMID(nf_fid,'z',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim z'
        return
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,z)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim z'
        return
      endif

      nf_status = nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
        return
      endif

      istatus = 1

      return
      end
C
C  Subroutine to read the file "LAPS fua/fsf files
C  - forecast model attributes" 
C
      subroutine read_lapsprd_attr(fullname,
     +     Dx, Dy, La1, Lo1, Latin1, Latin2, LoV, 
     +     grid_type, La2,Lo2, istatus)
C
      implicit none
      include 'netcdf.inc'
      integer nf_fid, nf_vid, nf_status
      integer Nx, Ny      !<-- not returned as these are x/y from dims call
      integer ln
      integer istatus
      real Dx, Dy, La1, La2, Lo1, Lo2, LoV
      real Latin1, Latin2
      character*132 grid_type_internal
      character*30 grid_type
      character*(*) fullname
C
C  Open netcdf File for reading
C
      istatus = 0
      call s_len(fullname,ln)
      nf_status = NF_OPEN(fullname,NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'NF_OPEN ',fullname(1:ln)
        return
      endif

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
C      La2          "last latitude"
C
      nf_status = NF_INQ_VARID(nf_fid,'La2',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var La1'
        return
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,La2)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var La2'
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
C      Lo2          "last longitude"
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
        print *,'in var Lo2'
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

C   Variables of type CHAR
C
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
      nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,grid_type_internal)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var grid_type'
        return 
      endif
      grid_type=grid_type_internal

      nf_status = nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
        return 
      endif

      istatus = 1
      return
      end
