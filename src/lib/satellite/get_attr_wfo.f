      subroutine get_attribute_wfo(filename,rlat00,rlon00,dx,
     &dy,nx,ny,istatus)
C
C  Open netcdf File for reading
C
      character filename*(*)
      character dummy*31
      integer   istatus
      integer   dim_id
      integer   nx,ny
      integer   rcode
      integer   ncid
      integer   lenf
      integer   nf_status
      real      rlat00
      real      rlon00
      real      dx,dy

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
C Get resolution-x
C
      nf_status = NF_INQ_ATTID(nf_fid,nf_attid,'dxKm',nf_attnum)
      if(nf_status.ne.NF_NOERR) then
         print*, NF_STRERROR(nf_status)
         print*, 'dxKm attribute id'
         istatus=-1
         goto 100
      endif

      nf_status = NF_GET_ATT_REAL(nf_fid,nf_attid,'dxKm',dx)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dxKm'
        istatus=-1
        goto 100
      endif
C
C Get resolution-y
C
      nf_status = NF_INQ_ATTID(nf_fid,nf_attid,'dyKm',nf_attnum)
      if(nf_status.ne.NF_NOERR) then
         print*, NF_STRERROR(nf_status)
         print*, 'dyKm attribute id'
         istatus=-1
         goto 100
      endif

      nf_status = NF_GET_ATT_REAL(nf_fid,nf_attid,'dyKm',dy)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dyKm'
        istatus=-1
        goto 100
      endif
c
c get x dimension
c
      dim_id = NCDID(nf_fid, 'x', rcode)
      if(rcode.ne.0)then
         write(6,*)'Error getting x id code - returning'
         istatus=-1
         goto 100
      endif

      call NCDINQ(nf_fid,dim_id,dummy,nx,RCODE)
      if(rcode.ne.0)then
         write(6,*)'Error getting x dimension - nx'
         istatus=-1
         goto 100
      endif
c
c get x dimension
c
      dim_id = NCDID(nf_fid, 'y', rcode)
      if(rcode.ne.0)then
         write(6,*)'Error getting y id code - returning'
         istatus=-1
         goto 100
      endif

      call NCDINQ(nf_fid,dim_id,dummy,ny,RCODE)
      if(rcode.ne.0)then
         write(6,*)'Error getting y dimension - ny'
         istatus=-1
         goto 100
      endif
C
      istatus = 0
100   return
      end
