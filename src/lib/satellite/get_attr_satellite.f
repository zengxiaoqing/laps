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
c     integer   ncid
      integer   lenf
      integer   nf_status
      real      rlat00
      real      rlon00
      real      dx,dy

      include 'netcdf.inc'

      lenf=index(filename,' ')-1
      print*,'opening ',filename(1:lenf)
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
      character c_afwa_fname*11

      cs1='i'
      if(chtype.eq.'vis')cs1='v'

      cs2=csatid(1:2)
      cs3=csatid(5:6)

      c_afwa_fname='u'//cs2//cs3//cs1//'1_'//chtype

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
