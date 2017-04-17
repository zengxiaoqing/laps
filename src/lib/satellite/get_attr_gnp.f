      subroutine get_attribute_gnp(nf_fid,centralLat,centralLon,
     &stdlat,stdlon,latNxNy,lonNxNy,latdxdy,londxdy,dx,dy,nx,ny,
     &istatus)
C
C  Open netcdf File for reading 
C
      character dummy*31
      integer   istatus
      integer   dim_id
      integer   nx,ny
      integer   rcode
      integer   nf_fid
      integer   lenf
      integer   nf_status
      integer   nf_attid
      integer   nf_attnum
      real      centralLat
      real      centralLon
      real      stdlat
      real      stdlon
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
C Get variable ID of projection parameters
      nf_status = NF_INQ_VARID(nf_fid,'lambert_projection',nf_projid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var lambert_projection'
      endif

C Get attributes, both global and from 'lambert_projection'      

C     
C Get centralLat
C
      nf_attid=nf_global
      nf_status = NF_INQ_ATTID(nf_fid,nf_attid,'tile_center_latitude'
     +                        ,nf_attnum)
      if(nf_status.ne.NF_NOERR) then
         print*, NF_STRERROR(nf_status)
         print*, 'tile_center_latitude attribute id'
         istatus=-1
         goto 100
      endif

      nf_status=NF_GET_ATT_REAL(nf_fid,nf_attid,'tile_center_latitude'
     +                         ,centralLat)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'tile_center_latitude'
         istatus=-1
         goto 100
      endif
C
C Get centralLon
C
      nf_attid=nf_global
      nf_status = NF_INQ_ATTID(nf_fid,nf_attid,'tile_center_longitude'
     +                        ,nf_attnum)
      if(nf_status.ne.NF_NOERR) then
         print*, NF_STRERROR(nf_status)
         print*, 'tile_center_longitude attribute id'
         istatus=-1
         goto 100
      endif

      nf_status=NF_GET_ATT_REAL(nf_fid,nf_attid,'tile_center_longitude'
     +                         ,centralLon)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'tile_center_longitude'
        istatus=-1
        goto 100
      endif
C
C Get standard_parallel
C
      nf_attid=nf_projid
      nf_status = NF_INQ_ATTID(nf_fid,nf_attid,'standard_parallel'
     +                        ,nf_attnum)
      if(nf_status.ne.NF_NOERR)	then
         print*, NF_STRERROR(nf_status)
         print*, 'standard_parallel attribute id'
         istatus=-1
         goto 100
      endif

      nf_status = NF_GET_ATT_REAL(nf_fid,nf_attid,'standard_parallel'
     +                           ,stdlat)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'standard_parallel'
         istatus=-1
         goto 100
      endif
C
C Get stdlon
C
      nf_attid=nf_projid
      nf_status = NF_INQ_ATTID(nf_fid,nf_attid
     +                       ,'longitude_of_central_meridian',nf_attnum)
      if(nf_status.ne.NF_NOERR) then
         print*, NF_STRERROR(nf_status)
         print*, 'longitude_of_central_meridian attribute id'
         istatus=-1
         goto 100
      endif

      nf_status = NF_GET_ATT_REAL(nf_fid,nf_attid
     +                          ,'longitude_of_central_meridian',stdlon)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'longitude_of_central_meridian'
        istatus=-1
        goto 100
      endif
C
C Get latNxNy
C
      nf_attid=nf_global
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
      nf_attid=nf_global
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
      nf_attid=nf_global
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
      nf_attid=nf_global
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
      nf_attid=nf_global
      nf_status = NF_INQ_ATTID(nf_fid,nf_attid,'pixel_x_size',nf_attnum)
      if(nf_status.ne.NF_NOERR) then
         print*, NF_STRERROR(nf_status)
         print*, 'pixel_x_size attribute id'
         istatus=-1
         goto 100
      endif

      nf_status = NF_GET_ATT_REAL(nf_fid,nf_attid,'pixel_x_size',dx)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'pixel_x_size'
        istatus=-1
        goto 100
      endif
C
C Get resolution-y
C
      nf_attid=nf_global
      nf_status = NF_INQ_ATTID(nf_fid,nf_attid,'pixel_y_size',nf_attnum)
      if(nf_status.ne.NF_NOERR) then
         print*, NF_STRERROR(nf_status)
         print*, 'pixel_y_size attribute id'
         istatus=-1
         goto 100
      endif

      nf_status = NF_GET_ATT_REAL(nf_fid,nf_attid,'pixel_y_size',dy)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'pixel_y_size'
        istatus=-1
        goto 100
      endif
c
c get x dimension
c
      nf_attid=nf_global
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
c get y dimension
c
      nf_attid=nf_global
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
