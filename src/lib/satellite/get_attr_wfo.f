      subroutine get_attribute_wfo(filename,rlat00,rlon00,resx,
     &resy,istatus)
C
C  Open netcdf File for reading
C
      character filename*(*)
      integer   istatus
      real      rlat00
      real      rlon00
      real      resx,resy

      include 'netcdf.inc'

      lenf=index(filename,' ')-1
      nf_status = NF_OPEN(filename,NF_NOWRITE,nf_fid)

      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'NF_OPEN ',filename(1:lenf)
        istatus=-1
        goto 100
      endif
C
C Get lat00
C
      nf_status = NF_INQ_ATTID(nf_fid,nf_attid,'lat00',nf_attnum)
      if(nf_status.ne.NF_NOERR)	then
         print*, NF_STRERROR(nf_status)
         print*, 'lat00 attribute id'
         istatus=-1
         goto 100
      endif

      nf_status = NF_GET_ATT_REAL(nf_fid,nf_attid,'lat00',rlat00)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'lat00'
         istatus=-1
         goto 100
      endif
C
C Get lon00
C
      nf_status = NF_INQ_ATTID(nf_fid,nf_attid,'lon00',nf_attnum)
      if(nf_status.ne.NF_NOERR) then
         print*, NF_STRERROR(nf_status)
         print*, 'lon00 attribute id'
         istatus=-1
         goto 100
      endif

      nf_status = NF_GET_ATT_REAL(nf_fid,nf_attid,'lon00',rlon00)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'lon00'
        istatus=-1
        goto 100
      endif
C
C Get resx
C
      nf_status = NF_INQ_ATTID(nf_fid,nf_attid,'dxKm',nf_attnum)
      if(nf_status.ne.NF_NOERR) then
         print*, NF_STRERROR(nf_status)
         print*, 'dxKm attribute id'
         istatus=-1
         goto 100
      endif

      nf_status = NF_GET_ATT_REAL(nf_fid,nf_attid,'dxKm',resx)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dxKm'
        istatus=-1
        goto 100
      endif
C
C Get resx
C
      nf_status = NF_INQ_ATTID(nf_fid,nf_attid,'dyKm',nf_attnum)
      if(nf_status.ne.NF_NOERR) then
         print*, NF_STRERROR(nf_status)
         print*, 'dyKm attribute id'
         istatus=-1
         goto 100
      endif

      nf_status = NF_GET_ATT_REAL(nf_fid,nf_attid,'dyKm',resy)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dyKm'
        istatus=-1
        goto 100
      endif
C
      istatus = 0
100   return
      end
