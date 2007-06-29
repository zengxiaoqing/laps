      subroutine get_bgdata_data
     +                   (i4time_sys,ilaps_cycle_time,NX_L,NY_L
     +                   ,i4time_earliest,i4time_latest
     +                   ,filename
     +                   ,lun_out
     +                   ,istatus)

      include 'netcdf.inc'

      character*(*) filename

      integer numIsoLevel, record, x, y,nf_fid, nf_vid, nf_status
C
C  Open netcdf File for reading
C
      nf_status=NF_OPEN(filename,NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),filename
        istatus=0
        return
      endif
C
C  Fill all dimension values
C
C
C Get size of numIsoLevel
C
      nf_status=NF_INQ_DIMID(nf_fid,'numIsoLevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim numIsoLevel'
      endif
      nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,numIsoLevel)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim numIsoLevel'
      endif
C
C Get size of record
C
      nf_status=NF_INQ_DIMID(nf_fid,'record',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim record'
      endif
      nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,record)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim record'
      endif
C
C Get size of x
C
      nf_status=NF_INQ_DIMID(nf_fid,'x',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim x'
      endif
      nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,x)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim x'
      endif
C
C Get size of y
C
      nf_status=NF_INQ_DIMID(nf_fid,'y',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim y'
      endif
      nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,y)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim y'
      endif

      nf_status=nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
      endif

      return
      end
C
C
C
      subroutine read_FMI_netcdf(filename, numIsoLevel, record, x, y, 
     +     isoLevel, DPT_2M, GH, GH_S, MSL, P, PVV, SH, T, T_2M, T_S,
     +     uW, uW_10M, vW, vW_10M, model, origin, reftime, valtime,
     +     istatus)
C
      include 'netcdf.inc'
      integer numIsoLevel, record, x, y,nf_fid, nf_vid, nf_status
      integer isoLevel(numIsoLevel)
      real DPT_2M( x,  y, record), GH( x,  y,  numIsoLevel, record),
     +     GH_S( x,  y, record), MSL( x,  y, record), P( x,  y,
     +     record), PVV( x,  y,  numIsoLevel, record), SH( x,  y, 
     +     numIsoLevel, record), T( x,  y,  numIsoLevel, record),
     +     T_2M( x,  y, record), uW( x,  y,  numIsoLevel, record),
     +     uW_10M( x,  y, record), vW( x,  y,  numIsoLevel, record),
     +     vW_10M( x,  y, record), T_S(x,  y,  record),
     +     SST(x,  y,  record), SKT(x,  y,  record), LSM(x, y, record)
      double precision reftime(record), valtime(record)
      character*132 origin
      character*132 model
      character*255 filename

      istatus =0
C
C  Open netcdf File for reading
C
      nf_status=NF_OPEN(filename,NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),TRIM(filename)
        istatus=1
        return
      endif

C   Variables of type REAL
C
C     Variable        NETCDF Long Name
C     DPT_2M        "2 metre dewpoint temperature"
C
      nf_status=NF_INQ_VARID(nf_fid,'DPT_2M',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for DPT_2M'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,DPT_2M)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for DPT_2M'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     GH            "geopotential"
C
      nf_status=NF_INQ_VARID(nf_fid,'GH',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for GH'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,GH)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for GH'
       endif
      endif
c convert to geometric height
c
      GH=GH/9.81
C
C     Variable        NETCDF Long Name
C     GH_S          "geopotential at surface"
C
      nf_status=NF_INQ_VARID(nf_fid,'GH_S',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for GH_S'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,GH_S)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for GH_S'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     MSL           "pressure at mean sea level"
C
      nf_status=NF_INQ_VARID(nf_fid,'MSL',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for MSL'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,MSL)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for MSL'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     P             "pressure at surface"
C
      nf_status=NF_INQ_VARID(nf_fid,'P',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for P'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,P)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for P'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     PVV           "Vertical velocity"
C
      nf_status=NF_INQ_VARID(nf_fid,'PVV',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for PVV'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,PVV)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for PVV'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     SH            "Specific humidity"
C
      nf_status=NF_INQ_VARID(nf_fid,'SH',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for SH'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,SH)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for SH'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     T             "temperature"
C
      nf_status=NF_INQ_VARID(nf_fid,'T',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for T'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,T)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for T'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     T_2M          "2 metre temperature"
C
      nf_status=NF_INQ_VARID(nf_fid,'T_2M',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for T_2M'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,T_2M)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for T_2M'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     LSM          "land/sea Mask"
C
      nf_status=NF_INQ_VARID(nf_fid,'LSM',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for LSM'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,LSM)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for LSM'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     SST          "Sea Surface Temperature"
C
      nf_status=NF_INQ_VARID(nf_fid,'SST',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for SST'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,SST)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for SST'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     SKT          "Skin temperature"
C
      nf_status=NF_INQ_VARID(nf_fid,'SKT',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for SKT'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,SKT)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for SKT'
       endif
      endif
C
C Prepare Surface Temperature array for return. SST for water
C points and skin temps for land points.
C
      do j=1,y
       do i=1,x
        if(LSM(i,j,1).eq.1)then
         T_S(i,j,1)=SST(i,j,1)
        else
         T_S(i,j,1)=SKT(i,j,1)
        endif
       enddo
      enddo
C
C     Variable        NETCDF Long Name
C     uW            "eastward wind"
C
      nf_status=NF_INQ_VARID(nf_fid,'uW',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for uW'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,uW)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for uW'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     uW_10M        "eastward wind"
C
      nf_status=NF_INQ_VARID(nf_fid,'uW_10M',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for uW_10M'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,uW_10M)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for uW_10M'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     vW            "northward wind"
C
      nf_status=NF_INQ_VARID(nf_fid,'vW',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for vW'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,vW)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for vW'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     vW_10M        "northward wind"
C
      nf_status=NF_INQ_VARID(nf_fid,'vW_10M',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for vW_10M'
      else
       nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,vW_10M)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for vW_10M'
       endif
      endif

C   Variables of type INT
C
C
C     Variable        NETCDF Long Name
C     isoLevel      "isobaric levels"
C
      nf_status=NF_INQ_VARID(nf_fid,'isoLevel',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for isoLevel'
      else
       nf_status=NF_GET_VAR_INT(nf_fid,nf_vid,isoLevel)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for isoLevel'
       endif
      endif

C   Variables of type DOUBLE
C
C
C     Variable        NETCDF Long Name
C     reftime       "reference time"
C
      nf_status=NF_INQ_VARID(nf_fid,'reftime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for reftime'
      else
       nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,reftime)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for reftime'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     valtime       "valid time"
C
      nf_status=NF_INQ_VARID(nf_fid,'valtime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for valtime'
      else
       nf_status=NF_GET_VAR_DOUBLE(nf_fid,nf_vid,valtime)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for valtime'
       endif
      endif


C   Variables of type CHAR
C
C
C     Variable        NETCDF Long Name
C     model         "model name"
C
      nf_status=NF_INQ_VARID(nf_fid,'model',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for model'
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,model)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for model'
       endif
      endif
C
C     Variable        NETCDF Long Name
C     origin        "origin name"
C
      nf_status=NF_INQ_VARID(nf_fid,'origin',nf_vid)
      if(nf_status.ne.NF_NOERR) then
       print *, NF_STRERROR(nf_status),' for origin'
      else
       nf_status=NF_GET_VAR_TEXT(nf_fid,nf_vid,origin)
       if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status),' for origin'
       endif
      endif

      nf_status=nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
      endif

      return
      end
