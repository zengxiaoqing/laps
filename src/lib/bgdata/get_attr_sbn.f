      subroutine get_attribute_sbn(cdfname,centralLat,centralLon,
     &rlat00,rlon00,latNxNy,lonNxNy,latdxdy,londxdy,dx,dy,nx,ny,
     &rotation,projname,istatus)
C
C  Open netcdf File for reading
C
      character cdfname*200
      character dummy*31
      character projname*30
      integer   istatus
      integer   dim_id
      integer   nx,ny
      integer   rcode
      integer   nf_fid
      integer   lenf,lenp
      integer   nf_status
      real      centralLat
      real      centralLon
      real      rlat00
      real      rlon00
      real      dx,dy
      real      latNxNy
      real      lonNxNy
      real      latdxdy
      real      londxdy
      real      rotation

      include 'netcdf.inc'

      print*,'In get_attribute_sbn'
      istatus = -1

      call s_len(cdfname,lenp)
      nf_status = NF_OPEN(cdfname,NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'NF_OPEN ', cdfname(1:lenp)
        return
      endif
C
C Get centralLat
C
      nf_status = NF_INQ_ATTID(nf_fid,nf_attid,'centralLat',nf_attnum)
      if(nf_status.ne.NF_NOERR) then
         print*, NF_STRERROR(nf_status)
         print*, 'centralLat attribute id'
         goto 100
      endif

      nf_status=NF_GET_ATT_REAL(nf_fid,nf_attid,'centralLat',centralLat)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'centralLat'
         goto 100
      endif
C
C Get centralLon
C
      nf_status = NF_INQ_ATTID(nf_fid,nf_attid,'centralLon',nf_attnum)
      if(nf_status.ne.NF_NOERR) then
         print*, NF_STRERROR(nf_status)
         print*, 'centralLon attribute id'
         goto 100
      endif

      nf_status=NF_GET_ATT_REAL(nf_fid,nf_attid,'centralLon',centralLon)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'centralLon'
        goto 100
      endif
C
C Get lat00
C
      nf_status = NF_INQ_ATTID(nf_fid,nf_attid,'lat00',nf_attnum)
      if(nf_status.ne.NF_NOERR)	then
         print*, NF_STRERROR(nf_status)
         print*, 'lat00 attribute id'
         goto 100
      endif

      nf_status = NF_GET_ATT_REAL(nf_fid,nf_attid,'lat00',rlat00)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'lat00'
         goto 100
      endif
C
C Get lon00
C
      nf_status = NF_INQ_ATTID(nf_fid,nf_attid,'lon00',nf_attnum)
      if(nf_status.ne.NF_NOERR) then
         print*, NF_STRERROR(nf_status)
         print*, 'lon00 attribute id'
         goto 100
      endif

      nf_status = NF_GET_ATT_REAL(nf_fid,nf_attid,'lon00',rlon00)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'lon00'
        goto 100
      endif
C
C Get latNxNy
C
      nf_status = NF_INQ_ATTID(nf_fid,nf_attid,'latNxNy',nf_attnum)
      if(nf_status.ne.NF_NOERR) then
         print*, NF_STRERROR(nf_status)
         print*, 'latNxNy attribute id'
         goto 100
      endif

      nf_status = NF_GET_ATT_REAL(nf_fid,nf_attid,'latNxNy',latNxNy)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'latNxNy'
         goto 100
      endif
C
C Get lonNxNy
C
      nf_status = NF_INQ_ATTID(nf_fid,nf_attid,'lonNxNy',nf_attnum)
      if(nf_status.ne.NF_NOERR) then
         print*, NF_STRERROR(nf_status)
         print*, 'lonNxNy attribute id'
         goto 100
      endif

      nf_status = NF_GET_ATT_REAL(nf_fid,nf_attid,'lonNxNy',lonNxNy)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'lonNxNy'
        goto 100
      endif
C
C Get latDxDy
C
      nf_status = NF_INQ_ATTID(nf_fid,nf_attid,'latDxDy',nf_attnum)
      if(nf_status.ne.NF_NOERR) then
         print*, NF_STRERROR(nf_status)
         print*, 'latDxDy attribute id'
         goto 100
      endif

      nf_status = NF_GET_ATT_REAL(nf_fid,nf_attid,'latDxDy',latdxdy)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'latDxDy'
         goto 100
      endif
C
C Get lonDxDy
C
      nf_status = NF_INQ_ATTID(nf_fid,nf_attid,'lonDxDy',nf_attnum)
      if(nf_status.ne.NF_NOERR) then
         print*, NF_STRERROR(nf_status)
         print*, 'lonDxDy attribute id'
         goto 100
      endif

      nf_status = NF_GET_ATT_REAL(nf_fid,nf_attid,'lonDxDy',londxdy)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'lonDxDy'
        goto 100
      endif
C
C Get resolution-x
C
      nf_status = NF_INQ_ATTID(nf_fid,nf_attid,'dxKm',nf_attnum)
      if(nf_status.ne.NF_NOERR) then
         print*, NF_STRERROR(nf_status)
         print*, 'dxKm attribute id'
         goto 100
      endif

      nf_status = NF_GET_ATT_REAL(nf_fid,nf_attid,'dxKm',dx)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dxKm'
        goto 100
      endif
C
C Get resolution-y
C
      nf_status = NF_INQ_ATTID(nf_fid,nf_attid,'dyKm',nf_attnum)
      if(nf_status.ne.NF_NOERR) then
         print*, NF_STRERROR(nf_status)
         print*, 'dyKm attribute id'
         goto 100
      endif

      nf_status = NF_GET_ATT_REAL(nf_fid,nf_attid,'dyKm',dy)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dyKm'
        goto 100
      endif
c
c get x dimension
c
      dim_id = NCDID(nf_fid, 'x', rcode)
      if(rcode.ne.0)then
         write(6,*)'Error getting x id code - returning'
         goto 100
      endif

      call NCDINQ(nf_fid,dim_id,dummy,nx,RCODE)
      if(rcode.ne.0)then
         write(6,*)'Error getting x dimension - nx'
         goto 100
      endif
c
c get x dimension
c
      dim_id = NCDID(nf_fid, 'y', rcode)
      if(rcode.ne.0)then
         write(6,*)'Error getting y id code - returning'
         goto 100
      endif

      call NCDINQ(nf_fid,dim_id,dummy,ny,RCODE)
      if(rcode.ne.0)then
         write(6,*)'Error getting y dimension - ny'
         goto 100
      endif

C
C Get resolution-y
C
	nf_status = NF_INQ_ATTID(nf_fid,nf_attid,'rotation',nf_attnum)
	if(nf_status.ne.NF_NOERR) then
	 print*, NF_STRERROR(nf_status)
	 print*, 'rotation: attribute id'
	 goto 100
	endif

	nf_status = NF_GET_ATT_REAL(nf_fid,nf_attid,'rotation'
     .,rotation)
	if(nf_status.ne.NF_NOERR) then
	print *, NF_STRERROR(nf_status)
	print *,'rotation'
	goto 100
	endif
C
C Get projection name
C
        nf_status = NF_INQ_ATTID(nf_fid,nf_attid,'projName',nf_attnum)
        if(nf_status.ne.NF_NOERR) then
         print*, NF_STRERROR(nf_status)
         print*, 'projName: attribute id'
         goto 100
        endif

        nf_status = NF_GET_ATT_TEXT(nf_fid,nf_attid,'projName'
     .,projname)
        if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'projname'
        goto 100
        endif

C
      istatus=1
100   return
      end
