      subroutine readavnpublicdims(fname,x,y,numIsoLevel,record,
     +istatus)

      include 'netcdf.inc'
      integer       numIsoLevel
      integer       record, x, y
      integer       nf_fid, nf_vid, nf_status
      integer       istatus
      character*200 fname
C
      istatus=0
C
C  Open netcdf File for reading
C
      nf_status = NF_OPEN(fname,NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        call s_len(fname,len)
        print *, NF_STRERROR(nf_status)
        print *,'NF_OPEN: ',fname(1:len)
        return
      endif
C
C  Fill all dimension values
C
C
C Get size of numIsoLevel
C
      nf_status = NF_INQ_DIMID(nf_fid,'numIsoLevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        nf_status = NF_INQ_DIMID(nf_fid,'isoLevel',nf_vid)
        if(nf_status.ne.NF_NOERR) then
           print *, NF_STRERROR(nf_status)
           print *,'dim numIsoLevel'
           return
        endif
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,numIsoLevel)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim numIsoLevel'
        return
      endif
C
C Get size of record
C
      nf_status = NF_INQ_DIMID(nf_fid,'record',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim record'
        return
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,record)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim record'
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

      nf_status = nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
        return
      endif

      istatus=1

C     call main_sub(nf_fid, numIsoLevel, record, x, y)

      return
      end
C
C *********************************************************************
C  Subroutine to read the file "Global 26 layer spectral Aviation Model" 
C
      subroutine read_avn_netcdf(fname, numIsoLevel, record, x, y, 
     +     version, ACCS, GH, GH_S, P, PMSL, PW, RH, T, T_2M, T_S,
     +     uW, vW, wW, RH_2M, isoLevel, reftime, valtime, grid, model,
     +     nav, origin, istatus)
C
      include 'netcdf.inc'
      integer numIsoLevel, record, x, y,nf_fid, nf_vid, nf_status
      integer version
      real ACCS( x,  y, record), GH( x,  y,  numIsoLevel, record),
     +     GH_S( x,  y, record), P( x,  y, record), PMSL( x,  y,
     +     record), PW( x,  y, record), RH( x,  y,  numIsoLevel,
     +     record), T( x,  y,  numIsoLevel, record), T_S( x,  y,
     +     record), T_2M( x,  y, record), uW( x,  y,  numIsoLevel,
     +     record), vW( x,  y, numIsoLevel, record), RH_2M( x, y,
     +     record), wW( x,  y, numIsoLevel, record)
      double precision isoLevel(numIsoLevel), reftime(record),
     +     valtime(record)
      character*132 origin
      character*132 model
      character*132 nav
      character*132 grid
      character*255 fname
C
C
      istatus = 1
C
C  Open netcdf File for reading
C
      nf_status = NF_OPEN(fname,NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        call s_len(fname,len)
        print *, NF_STRERROR(nf_status)
        print *,'NF_OPEN: ',fname(1:len)
        return
      endif
C
C     Variables of type REAL
C
C     Variable        NETCDF Long Name
C      ACCS         "accumulated snow"
C
      nf_status = NF_INQ_VARID(nf_fid,'ACCS',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var ACCS'
c       return
      else
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,ACCS)
        if(nf_status.ne.NF_NOERR) then
           print *, NF_STRERROR(nf_status)
           print *,'in var ACCS'
           return
        endif
      endif
C
C     Variable        NETCDF Long Name
C      GH           "geopotential height"
C
      nf_status = NF_INQ_VARID(nf_fid,'GH',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var GH'
        return
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,GH)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var GH'
        return
      endif
C
C     Variable        NETCDF Long Name
C      GH_S         "geopotential height of surface"
C
      nf_status = NF_INQ_VARID(nf_fid,'GH_S',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var GH_S'
        return
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,GH_S)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var GH_S'
        return
      endif
C
C     Variable        NETCDF Long Name
C      P            "pressure at surface"
C
      nf_status = NF_INQ_VARID(nf_fid,'P',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var P'
        return
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,P)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var P'
        return
      endif
C
C     Variable        NETCDF Long Name
C      PMSL         "pressure at mean sea level"
C
      nf_status = NF_INQ_VARID(nf_fid,'PMSL',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var PMSL'
        return
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,PMSL)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var PMSL'
        return
      endif
C
C     Variable        NETCDF Long Name
C      PW           "precipitable water"
C
      nf_status = NF_INQ_VARID(nf_fid,'PW',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var PW'
        return
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,PW)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var PW'
        return
      endif
C
C     Variable        NETCDF Long Name
C      RH           "relative humidity"
C
      nf_status = NF_INQ_VARID(nf_fid,'RH',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var RH'
        return
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,RH)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var RH'
        return
      endif
C
C     Variable        NETCDF Long Name
C      T            "temperature"
C
      nf_status = NF_INQ_VARID(nf_fid,'T',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var T'
        return
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,T)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var T'
        return
      endif
C
C     Variable        NETCDF Long Name
C      T_2M          "temperature at 2 meters above surface"
C
      nf_status = NF_INQ_VARID(nf_fid,'T_2M',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var T_2M'
        return
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,T_2M)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var T_2M'
        return
      endif
C
C     Variable        NETCDF Long Name
C      T_S          "temperature at surface"
C
      nf_status = NF_INQ_VARID(nf_fid,'T_S',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var T_S'
        return
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,T_S)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var T_S'
        return
      endif
C
C     Variable        NETCDF Long Name
C      uW           "eastward wind"
C
      nf_status = NF_INQ_VARID(nf_fid,'uW',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var uW'
        return
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,uW)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var uW'
        return
      endif
C
C     Variable        NETCDF Long Name
C      vW           "northward wind"

      nf_status = NF_INQ_VARID(nf_fid,'vW',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var vW'
        return
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,vW)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var vW'
        return
      endif
C
C     Variable        NETCDF Long Name
C      RH_2M          "Relative Humidity at 2 meters above surface"
C
      nf_status = NF_INQ_VARID(nf_fid,'RH_2M',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var RH_2M'
        return
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,RH_2M)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var RH_2M'
        return
      endif
C
C     Variable        NETCDF Long Name
C      wW           "pressure vertical velocity"

      nf_status = NF_INQ_VARID(nf_fid,'PVV',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var PVV (pres vert velocity)'
        return
      endif
      nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,wW)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var wW'
        return
      endif
C
C   Variables of type INT
C
C
C     Variable        NETCDF Long Name
C      version      
C
      nf_status = NF_INQ_VARID(nf_fid,'version',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var version'
        return
      endif
      nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,version)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var version'
        return
      endif

C   Variables of type DOUBLE
C
C
C     Variable        NETCDF Long Name
C      isoLevel     "isobaric levels"
C
      nf_status = NF_INQ_VARID(nf_fid,'isoLevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var isoLevel'
        return
      endif
      nf_status = NF_GET_VAR_DOUBLE(nf_fid,nf_vid,isoLevel)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var isoLevel'
        return
      endif
C
C     Variable        NETCDF Long Name
C      reftime      "reference time"
C
      nf_status = NF_INQ_VARID(nf_fid,'reftime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var reftime'
        return
      endif
      nf_status = NF_GET_VAR_DOUBLE(nf_fid,nf_vid,reftime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var reftime'
        return
      endif
C
C     Variable        NETCDF Long Name
C      valtime      "valid time"
C
      nf_status = NF_INQ_VARID(nf_fid,'valtime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var valtime'
        return
      endif
      nf_status = NF_GET_VAR_DOUBLE(nf_fid,nf_vid,valtime)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var valtime'
        return
      endif


C   Variables of type CHAR
C
C
C     Variable        NETCDF Long Name
C      grid_type 
C
      nf_status = NF_INQ_VARID(nf_fid,'grid_type',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var grid'
        return
      endif
      nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,grid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var grid'
        return
      endif
C
C     Variable        NETCDF Long Name
C      model        
C
      nf_status = NF_INQ_VARID(nf_fid,'model',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var model'
        return
      endif
      nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,model)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var model'
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
