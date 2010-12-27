      subroutine get_attribute_vol(nf_fid,StationLatitude
     &  ,StationLongitude,StationElevationInMeters,station 
!    &  ,latNxNy,lonNxNy,latdxdy,londxdy,dx,dy,nx,ny
     &  ,istatus)
C
C  Obtain attributes, still should add radar name
C
      character dummy*31
      character station*31
      integer   istatus
      integer   dim_id
      integer   nx,ny
      integer   rcode
      integer   nf_fid
      integer   lenf
      integer   nf_status
      integer   nf_attid
      integer   nf_attnum
      real      StationLatitude
      real      StationLongitude
      real      rStationElevationInMeters
      real      rlon00
      real      dx,dy
      real      latNxNy
      real      lonNxNy
      real      latdxdy
      real      londxdy

      include 'netcdf.inc'

c     lenf=index(filename,' ')-1
c     print*,'opening ',filename(1:lenf)
c     nf_status = NF_OPEN(filename,NF_NOWRITE,nf_fid)

c     if(nf_status.ne.NF_NOERR) then
c       print *, NF_STRERROR(nf_status)
c       print *,'NF_OPEN ',filename(1:lenf)
c       istatus=-1
c       goto 100
c     endif
C
C Get StationLatitude
C
      nf_attid=0
      nf_attnum=0
      nf_status = NF_INQ_ATTID(nf_fid,nf_attid,'StationLatitude'
     .                        ,nf_attnum)
      if(nf_status.ne.NF_NOERR) then
         print*, NF_STRERROR(nf_status)
         print*, 'StationLatitude attribute id'
         istatus=-1
         goto 100
      endif

      nf_status=NF_GET_ATT_REAL(nf_fid,nf_attid,'StationLatitude'
     .                         ,StationLatitude)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'StationLatitude'
         istatus=-1
         goto 100
      endif
C
C Get StationLongitude
C
      nf_status = NF_INQ_ATTID(nf_fid,nf_attid,'StationLongitude'
     .                        ,nf_attnum)
      if(nf_status.ne.NF_NOERR) then
         print*, NF_STRERROR(nf_status)
         print*, 'StationLongitude attribute id'
         istatus=-1
         goto 100
      endif

      nf_status=NF_GET_ATT_REAL(nf_fid,nf_attid,'StationLongitude'
     .                         ,StationLongitude)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'StationLongitude'
        istatus=-1
        goto 100
      endif
C
C Get StationElevationInMeters
C
      nf_status = NF_INQ_ATTID(nf_fid,nf_attid
     .                        ,'StationElevationInMeters',nf_attnum)
      if(nf_status.ne.NF_NOERR)	then
         print*, NF_STRERROR(nf_status)
         print*, 'StationElevationInMeters attribute id'
         istatus=-1
         goto 100
      endif

      nf_status = NF_GET_ATT_REAL(nf_fid,nf_attid
     .                           ,'StationElevationInMeters'
     .                           ,StationElevationInMeters)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'StationElevationInMeters'
         istatus=-1
         goto 100
      endif

C
C Get Station     
C
      nf_status = NF_INQ_ATTID(nf_fid,nf_attid,'Station',nf_attnum)
      if(nf_status.ne.NF_NOERR) then
         print*, NF_STRERROR(nf_status)
         print*, 'Station attribute id'
         istatus=-1
         goto 100
      endif

      nf_status = NF_GET_ATT_TEXT(nf_fid,nf_attid,'Station',station)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'Station'
        istatus=-1
        goto 100
      endif

      return
C
C Get latNxNy
C
      nf_status = NF_INQ_ATTID(nf_fid,nf_attid,'latNxNy',nf_attnum)
      if(nf_status.ne.NF_NOERR) then
         print*, NF_STRERROR(nf_status)
         print*, 'latNxNy attribute id'
         istatus=-1
         goto 100
      endif

      nf_status = NF_GET_ATT_REAL(nf_fid,nf_attid,'latNxNy',latNxNy)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'latNxNy'
         istatus=-1
         goto 100
      endif
C
C Get lonNxNy
C
      nf_status = NF_INQ_ATTID(nf_fid,nf_attid,'lonNxNy',nf_attnum)
      if(nf_status.ne.NF_NOERR) then
         print*, NF_STRERROR(nf_status)
         print*, 'lonNxNy attribute id'
         istatus=-1
         goto 100
      endif

      nf_status = NF_GET_ATT_REAL(nf_fid,nf_attid,'lonNxNy',lonNxNy)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'lonNxNy'
        istatus=-1
        goto 100
      endif
C
C Get latDxDy
C
      nf_status = NF_INQ_ATTID(nf_fid,nf_attid,'latDxDy',nf_attnum)
      if(nf_status.ne.NF_NOERR) then
         print*, NF_STRERROR(nf_status)
         print*, 'latDxDy attribute id'
         istatus=-1
         goto 100
      endif

      nf_status = NF_GET_ATT_REAL(nf_fid,nf_attid,'latDxDy',latdxdy)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'latDxDy'
         istatus=-1
         goto 100
      endif
C
C Get lonDxDy
C
      nf_status = NF_INQ_ATTID(nf_fid,nf_attid,'lonDxDy',nf_attnum)
      if(nf_status.ne.NF_NOERR) then
         print*, NF_STRERROR(nf_status)
         print*, 'lonDxDy attribute id'
         istatus=-1
         goto 100
      endif

      nf_status = NF_GET_ATT_REAL(nf_fid,nf_attid,'lonDxDy',londxdy)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'lonDxDy'
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
      istatus=0
100   return
      end
