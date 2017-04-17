      subroutine get_attribute_wfo(nf_fid,centralLat,centralLon,
     &rlat00,rlon00,latNxNy,lonNxNy,latdxdy,londxdy,dx,dy,nx,ny,
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
      real      rlat00
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
C Get centralLat
C
      nf_attid=0
      nf_attnum=0
      nf_status = NF_INQ_ATTID(nf_fid,nf_attid,'centralLat',nf_attnum)
      if(nf_status.ne.NF_NOERR) then
         print*, NF_STRERROR(nf_status)
         print*, 'centralLat attribute id'
         istatus=-1
         goto 100
      endif

      nf_status=NF_GET_ATT_REAL(nf_fid,nf_attid,'centralLat',centralLat)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'centralLat'
         istatus=-1
         goto 100
      endif
C
C Get centralLon
C
      nf_status = NF_INQ_ATTID(nf_fid,nf_attid,'centralLon',nf_attnum)
      if(nf_status.ne.NF_NOERR) then
         print*, NF_STRERROR(nf_status)
         print*, 'centralLon attribute id'
         istatus=-1
         goto 100
      endif

      nf_status=NF_GET_ATT_REAL(nf_fid,nf_attid,'centralLon',centralLon)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'centralLon'
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
c get y dimension
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
c
c============================================================
c
      function c_afwa_fname(csatid,chtype)

      character csatid*(*)
      character chtype*(*)
      character cs1*1
      character cs2*2
      character cs3*2
      character c_afwa_fname*(*)

      cs1='i'
      if(chtype.eq.'vis')cs1='v'

      if(csatid.eq.'meteos')then
         c_afwa_fname=csatid//cs1//'1_'//chtype
      else
         cs2=csatid(1:2)
         cs3=csatid(5:6)
         c_afwa_fname='u'//cs2//cs3//cs1//'1_'//chtype
      endif

      return 
      end

c
c ===========================================================
c
      function nw_vis_line_gwc(chtype,decimation,bescnfc,fsci)

      integer   nw_vis_line_gwc
      integer   fsci
      integer   bescnfc
      integer   decimation
      character chtype*3

      if(chtype.eq.'vis')then
         nw_vis_line_gwc=(bescnfc*decimation+fsci)
      else
         nw_vis_line_gwc=(bescnfc*decimation+fsci)*4 
      endif

      return
      end
c
c ===========================================================
c
      function nw_vis_pix_gwc(chtype,decimation,bepixfc,goalpha)

      integer    nw_vis_pix_gwc
      integer    bepixfc
      integer    decimation
      real       goalpha
      character  chtype*3

      if(chtype.eq.'vis')then
         nw_vis_pix_gwc=(bepixfc*decimation/4+int(goalpha))*8
      else
         nw_vis_pix_gwc=(bepixfc*decimation+int(goalpha))*8
      endif 

      return
      end
