      subroutine get_ruc2_dims(filename,cmodel,NX,NY,NZ
     &,StdLat1,StdLat2,Lon0,La1,Lo1,La2,Lo2,istatus)
c
c updated to get all nav info - JS 4-01
c
      implicit none
      include 'netcdf.inc'

      real    StdLat1,StdLat2
      real    Lon0
      real    La1,Lo1
      real    La2,Lo2

      integer  NX, NY, NZ, nf_fid, nf_vid, nf_status
      character filename*200
      character cmodel*(*)
      integer istatus

      istatus = 0
C
C  Open netcdf File for reading
C
      nf_status = NF_OPEN(filename,NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'NF_OPEN ', filename
        return
      endif
C
C  Fill all dimension values
C
C
C Get size of record
C
c      nf_status = NF_INQ_DIMID(nf_fid,'record',nf_vid)
c      if(nf_status.ne.NF_NOERR) then
c        print *, NF_STRERROR(nf_status)
c        print *,'dim record'
c      endif
c      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,record)
c      if(nf_status.ne.NF_NOERR) then
c        print *, NF_STRERROR(nf_status)
c        print *,'dim record'
c      endif
C
C Get size of x
C
      nf_status = NF_INQ_DIMID(nf_fid,'x',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim x'
        return
      endif

      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,NX)
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

      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,NY)
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

      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,NZ)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim z'
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
      if(cmodel.eq.'RUC40_NATIVE')then
        nf_status = NF_INQ_VARID(nf_fid,'Latin',nf_vid)
      else
        nf_status = NF_INQ_VARID(nf_fid,'Latin1',nf_vid)
      endif
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var Latin'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,StdLat1)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var Latin - StdLat'
      endif

c not in netcdf file
      StdLat2=StdLat1
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
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,Lon0)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var Lov - Lon0'
      endif

c these not in netcdf file!
      La2=55.4818
      Lo2=-57.3794

      nf_status = NF_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close ruc2'
        return
      endif

      istatus=1 
      return
      end
C
C
C
C  Subroutine to read the file "Hybrid-B 40km Rapid Update Cycle" 
C
      subroutine read_ruc2_hybb(filename, NX, NY, NZ, mmsp
     +     ,hgt, p, qv, u, v, vpt, w,istatus)
      implicit none
      include 'netcdf.inc'
C     integer NX, NY, NZ, nf_fid, nf_vid, nf_status,istatus
      integer NX, NY, NZ, nf_fid, nf_status,istatus
      character*200 filename
      character*4   cvar
      integer nxny,nxnynz
      real mmsp(nx,ny), hgt( NX,  NY,  NZ), 
     +     p( NX,  NY,  NZ), qv( NX,  NY,  NZ), 
     +     u( NX,  NY,  NZ), v( NX,  NY,  NZ), 
     +     vpt( NX,  NY,  NZ), w( NX,  NY,  NZ)

      istatus = 1
      nxny=nx*ny
      nxnynz=nx*ny*nz

C
C  Open netcdf File for reading
C
      nf_status = NF_OPEN(filename,NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'NF_OPEN ', filename
        return
      endif

C
C     Variable        NETCDF Long Name
C      MMSP         "MAPS mean sea level pressure" 
C
      cvar='MMSP'
      call read_netcdf_real(nf_fid,cvar,nxny,mmsp,0,0,nf_status)

C
C     Variable        NETCDF Long Name
C      hgt          "geopotential height" 
C
      cvar='hgt '
      call read_netcdf_real(nf_fid,cvar,nxnynz,hgt,0,0,nf_status)

      if(nf_status.lt.-0.5*nxnynz) then
         print*, 'A substantial portion of the height field is missing'
         print*, 'ABORTING file processing for ruc2 file ',filename
         return
      endif


C
C     Variable        NETCDF Long Name
C      p            "pressure" 
C
      cvar='p   '
      call read_netcdf_real(nf_fid,cvar,nxnynz,p,0,0,nf_status)

      if(nf_status.lt.-0.5*nxnynz) then
      print*, 'A substantial portion of the pressure field is missing'
         print*, 'ABORTING file processing for ruc2 file ',filename
         return
      endif
C
C     Variable        NETCDF Long Name
C      qv           "water vapor mixing ratio" 
C
      cvar='qv  '
      call read_netcdf_real(nf_fid,cvar,nxnynz,qv,0,0,nf_status)

      if(nf_status.lt.-0.5*nxnynz) then
      print*,'A substantial portion of the water vapor field is missing'
         print*, 'ABORTING file processing for ruc2 file ',filename
         return
      endif

C
C     Variable        NETCDF Long Name
C      u            "u-component of wind" 
C
      cvar='u   '
      call read_netcdf_real(nf_fid,cvar,nxnynz,u,0,0,nf_status)

      if(nf_status.lt.-0.5*nxnynz) then
         print*, 'A substantial portion of the u-wind field is missing'
         print*, 'ABORTING file processing for ruc2 file ',filename
         return
      endif

C
C     Variable        NETCDF Long Name
C      v            "v-component of wind" 
C
      cvar='v   '
      call read_netcdf_real(nf_fid,cvar,nxnynz,v,0,0,nf_status)

      if(nf_status.lt.-0.5*nxnynz) then
         print*, 'A substantial portion of the v-wind field is missing'
         print*, 'ABORTING file processing for ruc2 file ',filename
         return
      endif

C
C     Variable        NETCDF Long Name
C      vpt          "virtual potential temperature" 
C
      cvar='vpt '
      call read_netcdf_real(nf_fid,cvar,nxnynz,vpt,0,0,nf_status)

      if(nf_status.lt.-0.5*nxnynz) then
         print*, 'A substantial portion of the vpt field is missing'
         print*, 'ABORTING file processing for ruc2 file ',filename
         return
      endif

C
C     Variable        NETCDF Long Name
C      w            "vertical velocity" 
C
      cvar='w   '
      call read_netcdf_real(nf_fid,cvar,nxnynz,w,0,0,nf_status)


      nf_status = nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
        return
      endif

      istatus = 0
      return
      end
