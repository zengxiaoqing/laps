      subroutine get_unidata_model_id(filename,cmodel,ivaltimes,ntbg
     &,istatus)

      implicit none
      include 'netcdf.inc'
      character*132 cmodel
      character*200 filename
      character*132 model
      integer ntbg,istatus
      integer ivaltimes(ntbg)
      integer nf_fid,nf_vid,nf_status
C
C  Open netcdf File for reading
C
      istatus = 1
      nf_status = NF_OPEN(filename,NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'NF_OPEN ', filename
        return
      endif

      nf_status = NF_INQ_VARID(nf_fid,'model',nf_vid)
      if(nf_status.ne.NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *,'in var model'
         return
      endif
      nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,model)
      if(nf_status.ne.NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *,'in NF_GET_VAR_ model '
         return
      endif
      nf_status=NF_INQ_VARID(nf_fid,'valtimeMINUSreftime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *,'in NF_GET_VAR_ model '
         return
      endif
      nf_status=NF_GET_VARA_INT(nf_fid,nf_vid,1,ntbg,ivaltimes)
      if(nf_status.ne.NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *,'in NF_GET_VAR_ model '
         return
      endif

      nf_status = nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
        return
      endif

      if(model(1:3).ne.cmodel(1:3))then
         print*,'Mismatch between model and cmodel'
         print*,model,cmodel
      endif

      istatus = 0

      return
      end

      subroutine get_unidata_dims(cdfname,cmodel
     +,nxbg,nybg,nzbg_ht,nzbg_tp,nzbg_sh,nzbg_uv
     +,nzbg_ww,n_valtimes,istatus)

      implicit none
      include 'netcdf.inc'
      integer slen, nf_status,nf_fid, i, istat
      integer nf_vid
      character*132 cmodel
      character*200 cdfname

      integer nxbg,nybg
      integer nzbg_ht
      integer nzbg_tp
      integer nzbg_sh
      integer nzbg_uv
      integer nzbg_ww
      integer n_valtimes 
      integer record
      integer nvars
      integer ivaltimes(100)
      integer istatus

      integer ncid,itype,ndims
      integer j,k,kk,lc,nclen,lenc
      integer dimlen
      character*16 cvars(10)

      integer dimids(10)
      integer idims(10,10)
      integer nattr
      integer nf_attid,nf_attnum
      character*13 fname9_to_wfo_fname13, fname13
C
      istatus = 0
      call s_len(cdfname,slen)
C
C Get size of n_valtimes
C

      call get_nvaltimes_unidata(cdfname,n_valtimes,ivaltimes,istatus)
      if(istatus.ne.1) then
         print *,'Error: get_nvaltimes '
         return
      endif
C
C
C  Open netcdf File for reading
C
      nf_status = NF_OPEN(cdfname,NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'NF_OPEN unidata'
      endif
C
C Get everything for each variable
C
      call s_len(cmodel,nclen)
 
      nvars = 6 
      IF ( (cmodel(1:nclen) .EQ. 'RUC_ISO') .OR. 
     +     (cmodel(1:nclen) .EQ. 'GFS_ISO') ) THEN
        cvars(1)='Z               '
        cvars(2)='RH              '
        cvars(3)='T               '
        cvars(4)='u               ' 
        cvars(5)='v               '
        cvars(6)='omega           '
      ELSEIF(cmodel(1:nclen) .EQ. 'RUC_HYB')  THEN
        cvars(1)='Z_hybr          '
        cvars(2)='hum_mix_hybr    '
        cvars(3)='vptmp_hybr      '
        cvars(4)='u_hybr          '
        cvars(5)='v_hybr          '
        cvars(6)='omega_hybr      '
      ELSE
        PRINT *, "Unsupported Unidata model: ", cmodel
        istatus = 0
        return
      ENDIF

      do i=1,nvars

         nf_status = NF_INQ_VARID(nf_fid, cvars(i),nf_vid)
         nf_status = NF_INQ_VAR(nf_fid,nf_vid,cvars(i)
     +,itype,ndims,dimids,nattr)
         print *, cvars(i),nf_vid,ndims
         do j=1,ndims
            nf_status = NF_INQ_DIMLEN(nf_fid,dimids(j),dimlen)
            idims(j,i)= dimlen
         enddo

         if(i.eq.1)then
            nxbg = idims(1,i)
            nybg = idims(2,i)
            nzbg_ht=idims(3,i)
         elseif(i.eq.3)then
            nzbg_tp=idims(3,i)
         elseif(i.eq.2)then
            nzbg_sh=idims(3,i)
         elseif(i.eq.4)then
            nzbg_uv=idims(3,i)
         elseif(i.eq.6)then
            nzbg_ww=idims(3,i)
         endif

      enddo

      nf_status = nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
        return
      endif

      istatus = 1
      return 
      end

      subroutine get_nvaltimes_unidata(cdfname,nvaltimes,ivaltimes,
     +    istatus)

      implicit none

      include 'netcdf.inc'

      integer       nf_fid,nf_status,nf_vid
      integer, intent(out)::       nvaltimes
      integer, intent(out)::    ivaltimes(100)
      integer       ireftimes(100)
      integer,allocatable ::       ilocaltimes(:)
      integer,intent(out)::   istatus
      integer t
      character*200,intent(in) :: cdfname

c     logical       l2

      istatus=0

c     l2=.false.
c     inquire(file=cdfname,opened=l2)
c     if(.not.l2)then

      nf_status = NF_OPEN(cdfname,NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *,'NF_OPEN: get_nvaltimes'
      endif

c     endif
         
c switched from n_valtimes to record to handle incomplete files
c LW 7-9-03
c     nf_status = NF_INQ_DIMID(nf_fid,'n_valtimes',nf_vid)
      nf_status = NF_INQ_DIMID(nf_fid,'record',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim n_valtimes'
        return
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,nvaltimes)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim n_valtimes'
        return
      endif
      nf_status=NF_INQ_VARID(nf_fid,'valtime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *,'in NF_GET_VAR_ model '
         return
      endif
      ALLOCATE(ilocaltimes(nvaltimes))
      nf_status=NF_GET_VARA_INT(nf_fid,nf_vid,1,nvaltimes,ilocaltimes)
      if(nf_status.ne.NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *,'in NF_GET_VAR_ model '
         return
      endif

      ! Get the reftimes
      nf_status=NF_INQ_VARID(nf_fid,'reftime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *,'in NF_GET_VAR_ model '
         return
      endif
      nf_status=NF_GET_VARA_INT(nf_fid,nf_vid,1,nvaltimes,ireftimes)
      if(nf_status.ne.NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *,'in NF_GET_VAR_ model '
         return
      endif

      do t = 1, nvaltimes

        ivaltimes(t) = (ilocaltimes(t) - ireftimes(t)) * 3600
      enddo
      DEALLOCATE(ilocaltimes)
      nf_status = nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *,'nf_close: get_nvaltimes'
         return
      endif

      istatus=1
      return
      end
C
C ------------------------------------------------------------
      subroutine read_unidata_iso(cdfname,af,cmodel,
     .nxbg,nybg,nzbght,nzbgtp,nzbgsh,nzbguv,nzbgww,
     .prbght,prbgsh,prbguv,prbgww,
     .ht,tp,sh,uw,vw,ww,
     .ht_sfc,pr_sfc,uw_sfc,vw_sfc,sh_sfc,tp_sfc,mslp,
     .ctype,istatus)
c
      implicit none
c
      include 'netcdf.inc'
      include 'bgdata.inc'

c     integer ncid, lenstr, ntp, nvdim, nvs, ndsize
      integer model_out
      integer ncid

c     model_out=1  => lga
c     model_out=2  => dprep

      integer ndims ,dimids(NF_MAX_VAR_DIMS)
      integer itype,nattr

      integer nxbg,nybg
      integer nzbght
      integer nzbgtp
      integer nzbgsh
      integer nzbguv
      integer nzbgww
      integer nzunidata
      integer ntbg
      integer rcode
      integer ivaltimes(100)
      integer ind2m, ind10m
      logical lcmpsfcq
c
      real, intent(out)  ::   mslp(nxbg,nybg)

c *** 3D Output arrays.
c
      real, intent(out)  :: prbght(nxbg,nybg,nzbght)
      real, intent(out)  :: prbgsh(nxbg,nybg,nzbgsh)
      real, intent(out)  :: prbguv(nxbg,nybg,nzbguv)
      real, intent(out)  :: prbgww(nxbg,nybg,nzbgww)
      real, intent(out)  ::     ht(nxbg,nybg,nzbght)
      real, intent(out)  ::     tp(nxbg,nybg,nzbgtp)
      real, intent(out)  ::     sh(nxbg,nybg,nzbgsh)
      real, intent(out)  ::     uw(nxbg,nybg,nzbguv)
      real, intent(out)  ::     vw(nxbg,nybg,nzbguv)
      real, intent(out)  ::     ww(nxbg,nybg,nzbgww)
      real, intent(out)  ::     ht_sfc(nxbg,nybg)
      real, intent(out)  ::     tp_sfc(nxbg,nybg)
      real, intent(out)  ::     sh_sfc(nxbg,nybg)
      real, intent(out)  ::     uw_sfc(nxbg,nybg)
      real, intent(out)  ::     vw_sfc(nxbg,nybg)
      real, intent(out)  ::     pr_sfc(nxbg,nybg)
 
c
      real              ::  prbg(nzbght)

      integer start(10),count(10)
 
      integer i,j,k,n,ip,jp,ii,jj,it,kk
      integer istatus,slen,lent
      integer ibdht,ibdtp,ibduv,ibdsh,ibdww
c
      character*9   fname,oldfname,model
      character*5   ctype
      character*4   af
      character*16  cvar
      character*2   gproj
      character*200 cdfname
      character*132 cmodel
c
      real   xe,mrsat
      real   make_ssh

      integer nf_vid,nn,nf_status
      real cp,rcp, factor
      parameter (cp=1004.,rcp=287./cp)
c
c_______________________________________________________________________________
c
      interface

        subroutine read_netcdf_real(nf_fid,fname,n1,f
     +,start,count,istatus)
          integer n1
          integer nf_fid
          integer istatus
          integer start(10),count(10)
          real    f(n1)
          character*(*) fname
        end subroutine
      end interface
c
c -------------------------------------------------------

      print*,'HERE: read_unidata_iso'

      istatus = 1

      call s_len(cdfname,slen)

      print*,'cdfname: ',cdfname(1:slen)

      call get_nvaltimes_unidata(cdfname,ntbg,ivaltimes,istatus)
c
      print*,'opening cdf file: ',cdfname(1:slen)

      rcode = NF_OPEN(cdfname,NF_NOWRITE,ncid)
      if(rcode.ne.NF_NOERR) then
         print *, NF_STRERROR(rcode)
         print *,'NF_OPEN ',cdfname(1:slen)
         return
      endif

      read(af,'(i4)') nn

      n=1
      do while(n.lt.ntbg.and.ivaltimes(n)/3600.ne. nn)
         n=n+1
      enddo
      if(ivaltimes(n)/3600.ne.nn) then

         print*,'ERROR: No record valid at requested time '
         print*,'ntbg/nn/af/n/ivaltimes(n) ',ntbg,' ',nn,' ',af,
     &' ',n,' ',ivaltimes(n)

         rcode= NF_CLOSE(ncid)
         if(rcode.ne.NF_NOERR) then
            print *, NF_STRERROR(rcode)
            print *,'in NF_GET_VAR: ',cmodel
            return
         endif

         goto 999

      else

         print*,'Found valid record at ivaltime'
         print*,'ntbg/nn/af/n/ivaltimes(n) ',ntbg,' ',nn,' ',af,
     &' ',n,' ',ivaltimes(n) 
         print*
      endif

      ! Get the pressure levels for this data
      nf_status = NF_INQ_VARID(ncid,'level',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in level '
        istatus = 0
        return
      endif
      nf_status = NF_GET_VAR_REAL(ncid,nf_vid,prbg)
     
      ! Get index for 2m and 10m winds
      ind2m = 1
      ind10m = 2
  
      ! Set some indices
      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=1
      count(3)=nzbght
      start(4)=n
      count(4)=1

      print*,'read ht'
      cvar='Z'
      call read_netcdf_real(ncid,cvar,nxbg*nybg*count(3),ht,start
     +     ,count,rcode)
      if(rcode.ne.NF_NOERR) then
         if(rcode.gt.-61)then
            print *, NF_STRERROR(rcode)
            print *,'in NF_GET_VAR (Z): ',cmodel
         else
            print *,'Missing HT data detected: return'
         endif
         print*
         return
      endif

c
c ****** Statements to fill tp.
c
      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=1
      count(3)=nzbgtp
      start(4)=n
      count(4)=1

      print*,'read tp'
      cvar='T'
      call read_netcdf_real(ncid,cvar,nxbg*nybg*count(3),tp,start
     +     ,count,rcode)

      if(rcode.ne.NF_NOERR) then
         if(rcode.gt.-61)then
            print *, NF_STRERROR(rcode)
            print *,'in NF_GET_VAR (t): ',cmodel
         else
            print *,'Missing T data detected: return'
         endif
         print*
         return
      endif

c
c ****** Statements to fill rh.                           
c
      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=1
      count(3)=nzbgsh
      start(4)=n
      count(4)=1

      print*,'read rh'
      cvar='RH'
      call read_netcdf_real(ncid,cvar,nxbg*nybg*count(3),sh,start
     +     ,count,rcode)
     
      if(rcode.ne.NF_NOERR) then
         if(rcode.gt.-61)then
            print *, NF_STRERROR(rcode)
            print *,'in NF_GET_VAR (rh): ',cmodel
         else
c          RH is missing above 100 mb in GFS
           IF (cmodel .EQ. 'GFS_ISO') THEN
            print *, "FILLING RH AT TOP LEVELS!"
             do j=1,nybg
             do i=1,nxbg
             do k = 1,nzbgsh
               if ((sh(i,j,k) .LT. 0.).OR. 
     +             (sh(i,j,k) .GT. 200.)) THEN
                 sh(i,j,k) = 1.0
               endif
             enddo
             enddo
             enddo
           else
             print *,'Missing RH data detected: return'
             return
           endif
         endif
      endif
c
c ****** Statements to fill uw. 
c
      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=1
      count(3)=nzbguv
      start(4)=n
      count(4)=1

      print*,'read uw'
      cvar='u'
      call read_netcdf_real(ncid,cvar,nxbg*nybg*count(3),uw,start
     +     ,count,rcode)

      if(rcode.ne.NF_NOERR) then
         if(rcode.gt.-61)then
            print *, NF_STRERROR(rcode)
            print *,'in NF_GET_VAR (uw): ',cmodel
         else
            print *,'Missing U data detected: return'
         endif
         print*
         return
      endif

c
c ****** Statements to fill vw.                           
c
      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=1
      count(3)=nzbguv
      start(4)=n
      count(4)=1

      print*,'read vw'
      cvar='v'
      call read_netcdf_real(ncid,cvar,nxbg*nybg*count(3),vw,start
     +     ,count,rcode)

      if(rcode.ne.NF_NOERR) then
         if(rcode.gt.-61)then
            print *, NF_STRERROR(rcode)
            print *,'in NF_GET_VAR (vw): ',cmodel
         else
            print *,'Missing V data detected: return'
         endif
         print*
         return
      endif
c
c ****** Statements to fill ww.
c
      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=1
      count(3)=nzbgww
      start(4)=n
      count(4)=1

      print*,'read ww'
      cvar='omega'
      call read_netcdf_real(ncid,cvar,nxbg*nybg*count(3),ww  
     +  ,start,count,rcode)

      if(rcode.ne.NF_NOERR) then
         if(rcode.gt.-61)then
            print *, NF_STRERROR(rcode)
            print *,'in NF_GET_VAR (ww): ',cmodel
         else
C           ww is missing above 100 mb in GFS
           IF (cmodel .EQ. 'GFS_ISO') THEN
             do j=1,nybg
             do i=1,nxbg
             do k = 1,nzbgww
               if ((ww(i,j,k) .LT. -1000.) .OR.
     +             (ww(i,j,k) .GT. 1000)) THEN
                 ww(i,j,k) = 0.
               endif
             enddo
             enddo
             enddo
           else
             print *,'Missing WW data detected: return'
             print *,'Missing ww data detected: continue without'
             print *,'Filling ww with 0.0'
             ww(:,:,:) = 0.0
           endif
         endif
      endif

C   Get 2m T and RH, 10m U and V
      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=ind2m
      count(3)=1
      start(4)=n
      count(4)=1
      print*,'read tp_sfc'
      lcmpsfcq=.true.
      cvar='T_fhg'
      call read_netcdf_real(ncid,cvar,nxbg*nybg,tp_sfc,start
     +     ,count,rcode)
      if(rcode.ne.NF_NOERR) then
         if(rcode.gt.-61)then
            print *, NF_STRERROR(rcode)
            print *,'in NF_GET_VAR (tp_sfc): ',cmodel
         endif
      endif
      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=1
      count(3)=1
      start(4)=n
      count(4)=1
      print*,'read ht_sfc'
      lcmpsfcq=.true.
      cvar='Z_sfc'
      call read_netcdf_real(ncid,cvar,nxbg*nybg,ht_sfc,start
     +     ,count,rcode)
      if(rcode.ne.NF_NOERR) then
         if(rcode.gt.-61)then
            print *, NF_STRERROR(rcode)
            print *,'in NF_GET_VAR (ht_sfc): ',cmodel
         endif
      endif

      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=ind2m
      count(3)=1
      start(4)=n
      count(4)=1
      print*,'read sh_sfc'
      lcmpsfcq=.true.
      cvar='RH_fhg'
      call read_netcdf_real(ncid,cvar,nxbg*nybg,sh_sfc,start
     +     ,count,rcode)
      if(rcode.ne.NF_NOERR) then
         if(rcode.gt.-61)then
            print *, NF_STRERROR(rcode)
            print *,'in NF_GET_VAR (sh_sfc): ',cmodel
         endif
      endif

      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=ind10m
      count(3)=1
      start(4)=n
      count(4)=1
      print*,'read uw_sfc'
      lcmpsfcq=.true.
      cvar='u_fhg'
      call read_netcdf_real(ncid,cvar,nxbg*nybg,uw_sfc,start
     +     ,count,rcode)
      if(rcode.ne.NF_NOERR) then
         if(rcode.gt.-61)then
            print *, NF_STRERROR(rcode)
            print *,'in NF_GET_VAR (uw_sfc): ',cmodel
         endif
      endif

      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=ind10m
      count(3)=1
      start(4)=n
      count(4)=1
      print*,'read vw_sfc'
      lcmpsfcq=.true.
      cvar='v_fhg'
      call read_netcdf_real(ncid,cvar,nxbg*nybg,vw_sfc,start
     +     ,count,rcode)
      if(rcode.ne.NF_NOERR) then
         if(rcode.gt.-61)then
            print *, NF_STRERROR(rcode)
            print *,'in NF_GET_VAR (vw_sfc): ',cmodel
         endif
      endif

c
c get sfc pressure field
c
      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=1
      count(3)=1
      start(4)=n
      count(4)=1

      print*,'read p'

      lcmpsfcq=.true.
      
      cvar='P_sfc'
      call read_netcdf_real(ncid,cvar,nxbg*nybg,pr_sfc,start
     +     ,count,rcode)

      if(rcode.ne.NF_NOERR) then
         if(rcode.gt.-61)then
            print *, NF_STRERROR(rcode)
            print *,'in NF_GET_VAR (p): ',cmodel
         endif
         print *,'Missing sfc p data detected'
         print*,' -> continue without; compute in sfcbkgd'
         lcmpsfcq=.false.
         print*
c        return
      endif
c
c get mslp (this field name differs from one model to the other)
c
      if(cmodel(1:3).eq.'RUC')THEN
         cvar='Pm_msl'
      else
         cvar='P_msl'
      endif
      call read_netcdf_real(ncid,cvar,nxbg*nybg,mslp
     +           ,start,count,rcode)
      if(rcode.ne.NF_NOERR) then
        if(rcode.gt.-61)then
          print *, NF_STRERROR(rcode)
          print *,'in NF_GET_VAR (emsp): ',cmodel
        else
          print*,'Error status returned from read_netcdf_real'
        endif
        print *,'Missing emsp data detected'
        print*
        return
      endif

c
c *** Close netcdf file.
c
      rcode= NF_CLOSE(ncid)
      if(rcode.ne.NF_NOERR) then
         print *, NF_STRERROR(rcode)
         print *,'in NF_GET_VAR: ',cmodel
         return
      endif
c
ccc      endif

c
c *** Fill ouput arrays.
c *** Convert rh to sh.
c
      print*,'load prbg arrays'
      do j=1,nybg
      do i=1,nxbg
         prbgsh(i,j,:)=prbg(:)
         prbght(i,j,:)=prbg(:)
         prbguv(i,j,:)=prbg(:)
         prbgww(i,j,:)=prbg(:)
      enddo
      enddo

c for laps-lgb

      call s_len(ctype,lent)

      print*,'ctype ',ctype(1:lent)

c since we may not have pr_sfc at this point lets not
c do this comp here and instead let sfcbkgd do it later.
      if(ctype(1:lent).eq.'lapsb' .and. lcmpsfcq)then

         print*,' computing sfc q '
         ibdtp=0
         ibduv=0
         do j=1,nybg
         do i=1,nxbg
            if(tp_sfc(i,j).lt.missingflag.and.
     .         tp_sfc(i,j).gt.150.)       then
c
c make sfc q from rh
c
              it=int(tp_sfc(i,j)*100)
              it=min(45000,max(15000,it))
              xe=esat(it)
              mrsat=0.00622*xe/(pr_sfc(i,j)*0.01-xe)
              sh_sfc(i,j)=sh_sfc(i,j)*mrsat
              sh_sfc(i,j)=sh_sfc(i,j)/(1.+sh_sfc(i,j))
            else
              ibdtp=ibdtp+1
            endif
            if(uw_sfc(i,j).gt.500..or.vw_sfc(i,j).gt.500.)then
              ibduv=ibduv+1
            endif
         enddo
         enddo

         print*,'done computing sfc q'

      endif


      if(ibdtp.gt.0.or.ibduv.gt.0)then
         print*,'Found bad sfc data (tp/uv) ',ibdtp,ibduv
         return
      endif


      ibdht=0
      ibdtp=0
      ibduv=0
      do i=1,nxbg
      do j=1,nybg

         do k=1,nzbght              
            if (ht(i,j,k) .ge. 99999.)then
                ibdht=ibdht+1
            endif
         enddo
         do k=1,nzbgtp
            if(tp(i,j,k).lt.100.or.tp(i,j,k).ge.missingflag)then
               ibdtp=ibdtp+1
            endif
         enddo
         do k=1,nzbguv
            if(abs(uw(i,j,k)).ge. 500. .or.
     .         abs(vw(i,j,k)).ge. 500.)then
                ibduv=ibduv+1
            endif
         enddo
      enddo
      enddo

      if(ibdht.gt.0.or.ibdtp.gt.0.or.ibduv.gt.0)then
         print*,'Found bad 3d data (ht/t/uv) ',ibdht,ibdtp,ibduv
         print*,'Return to read_bgdata'
         return
      endif

      ibdsh=0
      if(ctype(1:lent).eq.'lapsb')then
         print*,'computing 3d sh'
         print*,'nzbgsh/nzbgtp = ',nzbgsh,nzbgtp
         do j=1,nybg
         do i=1,nxbg
         do k=1,nzbgsh
           if (sh(i,j,k).lt.200)then
             sh(i,j,k)= make_ssh(prbgsh(i,j,k)
     .,tp(i,j,k)-273.15,sh(i,j,k)/100., -132.)*0.001  !kg/kg 

           else
            ibdsh=ibdsh+1
           endif
         enddo
         enddo
         enddo

      endif

      print*,'done computing 3d sh'

      if(ibdsh.gt.0)then
         print*,'Found bad rh 3d data',ibdsh
         print*,'Return to read_bgdata'
         return
      endif


c this for dprep ... ingnore for now!
c -----------------------------------
      if(.false. .and. model_out.eq.2) then
c
c Compute exner and convert temp to theta
c
         do k=1,nzbght
            do j=1,nybg
               do i=1,nxbg
                  if(tp(i,j,k).ne.missingflag) then
                     factor=(1000./prbght(i,j,k))**rcp
                     tp(i,j,k) = tp(i,j,k)*factor
                     prbght(i,j,k) = cp/factor
                  endif
               enddo
            enddo
         enddo

      endif


c     if(istatus_211 .eq. 0) then
c       print*, 'No valid data found for',fname, af
c       return
c     endif

c      istatus = 1

      if(0.eq.1) then
 900     print*,'ERROR: bad dimension specified in netcdf file'
         print*, (count(i),i=1,4)
         istatus=-1
      endif
 999  return
      end

c -----------------------------------------------------------

      subroutine get_unidata_grid(filename,cmodel,NX,NY
     &,StdLat1,StdLat2,Lov,La1,Lo1,La2,Lo2,Dx,Dy,gproj, istatus)
c
      implicit none
      include 'netcdf.inc'

      real    StdLat1,StdLat2
      real    Lov
      real    La1,Lo1
      real    La2,Lo2
      real    Dx, Dy
      character*2 gproj
      character*100 grid_type
      integer  NX, NY, nf_fid, nf_vid, nf_status
      character filename*200
      character cmodel*(*)
      character*2, dimname
      integer istatus

      istatus = 0
C
C  Open netcdf File for reading
C
      print *, "OPENING:  ", trim(filename)
      nf_status = NF_OPEN(filename,NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'NF_OPEN ', filename
        return
      endif
C
C     Variable        NETCDF Long Name
C     grid_type
C
      nf_status = NF_INQ_VARID(nf_fid,'grid_type',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var grid_type'
      endif
      nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,grid_type)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var grid_type'
      endif
                                                                            
      IF (grid_type(1:7) .EQ. "Lambert") THEN
        gproj = "LC"
      ELSEIF(grid_type(1:18) .EQ. "Latitude/Longitude") THEN
        gproj = "LL"
      ELSE
        print *, "Found unsupported grid_type:",trim(grid_type)
        print *, "In get_unidata_grid"
        istatus = 0
        return
      ENDIF

C     Nx
      IF (gproj .NE. "LL") THEN
        dimname = "Nx"
      ELSE
        dimname = "Ni"
      ENDIF
      nf_status = NF_INQ_VARID(nf_fid,dimname,nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var ',dimname
        istatus = 0
        return
      else  
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,nx)
      endif

C     Ny
      IF (gproj .NE. "LL") THEN
        dimname = "Ny"
      ELSE
        dimname = "Nj"
      ENDIF

      nf_status = NF_INQ_VARID(nf_fid,dimname,nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var ',dimname
        istatus = 0
        return
      else
        nf_status = NF_GET_VAR_INT(nf_fid,nf_vid,ny)
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
      IF (lo1 .GT. 180.) lo1 = lo1 - 360.
C
C     Variable        NETCDF Long Name
C      Lov          "orientation of grid"
C
      IF (gproj .NE. "LL") THEN
        nf_status = NF_INQ_VARID(nf_fid,'Lov',nf_vid)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var Lov'
        endif
          nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,lov)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var Lov - lov'
        endif
        if (lov .GT. 180.) lov = lov - 360.
C
Ci      Variable        NETCDF Long Name
C       Latin1
C
        nf_status = NF_INQ_VARID(nf_fid,'Latin1',nf_vid)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var Latin1'
        endif
          nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,StdLat1)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var Latin1 - StdLat1'
        endif
C
C       Variable        NETCDF Long Name
C         Latin2
C
        nf_status = NF_INQ_VARID(nf_fid,'Latin2',nf_vid)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var Latin2'
        endif
          nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,StdLat2)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var Latin2 - StdLat2'
        endif
      ENDIF
C
C     Variable        NETCDF Long Name
C     Dx 
C
      IF (gproj .NE. "LL") THEN
        dimname = "Dx"
      ELSE 
        dimname = "Di"
      ENDIF 
      nf_status = NF_INQ_VARID(nf_fid,dimname,nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var Dx'
      endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,dx)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var Dx'
      endif
C
C     Variable        NETCDF Long Name
C     Dy
C
      IF (gproj .NE. "LL") THEN
        dimname = "Dy"
      ELSE
        dimname = "Dj"
      ENDIF
      nf_status = NF_INQ_VARID(nf_fid,dimname,nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var Dy'
      endif
       nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,dy)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in var Dy'
      endif

      ! Set or compute corner points

      IF (cmodel(1:3) .EQ. "RUC") THEN
        La2=55.481
        Lo2=-57.381
      ELSEIF(gproj .EQ. "LL") THEN
C
C       Variable        NETCDF Long Name
C       La2          "last latitude"
C
        nf_status = NF_INQ_VARID(nf_fid,'La2',nf_vid)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var La2'
        endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,La2)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var La2'
        endif
C
C       Variable        NETCDF Long Name
C       Lo2          "last longitude"
C
        nf_status = NF_INQ_VARID(nf_fid,'Lo2',nf_vid)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var Lo1'
        endif
        nf_status = NF_GET_VAR_REAL(nf_fid,nf_vid,Lo2)
        if(nf_status.ne.NF_NOERR) then
          print *, NF_STRERROR(nf_status)
          print *,'in var Lo1'
        endif
        IF (lo2 .GT. 180.) lo2 = lo2 - 360.

      ELSE
        print *, "Need to add computation of la2,lo2"
        print *, " In get_unidata_grid"
        istatus = 0
        return
      ENDIF

      nf_status = NF_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        return
      endif

      istatus=1 
      return
      end
C
C ------------------------------------------------------------
      subroutine read_unidata_ruc_hyb(cdfname,af,cmodel,
     .nxbg,nybg,nzbght,nzbgtp,nzbgsh,nzbguv,nzbgww,
     .prbght,prbgsh,prbguv,prbgww,
     .ht,tp,sh,uw,vw,ww,
     .ht_sfc,pr_sfc,uw_sfc,vw_sfc,sh_sfc,tp_sfc,mslp,
     .ctype,istatus)
c
      implicit none
c
      include 'netcdf.inc'
      include 'bgdata.inc'

c     integer ncid, lenstr, ntp, nvdim, nvs, ndsize
      integer model_out
      integer ncid

c     model_out=1  => lga
c     model_out=2  => dprep

      integer ndims ,dimids(NF_MAX_VAR_DIMS)
      integer itype,nattr

      integer nxbg,nybg
      integer nzbght
      integer nzbgtp
      integer nzbgsh
      integer nzbguv
      integer nzbgww
      integer nzunidata
      integer ntbg
      integer rcode
      integer ivaltimes(100)
      integer ind2m, ind10m
      logical lcmpsfcq
c
      real, intent(out)  ::   mslp(nxbg,nybg)

c *** 3D Output arrays.
c
      real, intent(out)  :: prbght(nxbg,nybg,nzbght)
      real, intent(out)  :: prbgsh(nxbg,nybg,nzbgsh)
      real, intent(out)  :: prbguv(nxbg,nybg,nzbguv)
      real, intent(out)  :: prbgww(nxbg,nybg,nzbgww)
      real, intent(out)  ::     ht(nxbg,nybg,nzbght)
      real, intent(out)  ::     tp(nxbg,nybg,nzbgtp)
      real, intent(out)  ::     sh(nxbg,nybg,nzbgsh)
      real, intent(out)  ::     uw(nxbg,nybg,nzbguv)
      real, intent(out)  ::     vw(nxbg,nybg,nzbguv)
      real, intent(out)  ::     ww(nxbg,nybg,nzbgww)
      real, intent(out)  ::     ht_sfc(nxbg,nybg)
      real, intent(out)  ::     tp_sfc(nxbg,nybg)
      real, intent(out)  ::     sh_sfc(nxbg,nybg)
      real, intent(out)  ::     uw_sfc(nxbg,nybg)
      real, intent(out)  ::     vw_sfc(nxbg,nybg)
      real, intent(out)  ::     pr_sfc(nxbg,nybg)
 
c
      integer start(10),count(10)
 
      integer i,j,k,n,ip,jp,ii,jj,it,kk
      integer istatus,slen,lent
      integer ibdht,ibdtp,ibduv,ibdsh,ibdww
c
      character*9   fname,oldfname,model
      character*5   ctype
      character*4   af
      character*16  cvar
      character*2   gproj
      character*200 cdfname
      character*132 cmodel
c
      real tv
      integer nf_vid,nn,nf_status,ic,jc
      real cp,rcp, factor
      parameter (cp=1004.,rcp=287./cp)
c
c_______________________________________________________________________________
c
      interface

        subroutine read_netcdf_real(nf_fid,fname,n1,f
     +,start,count,istatus)
          integer n1
          integer nf_fid
          integer istatus
          integer start(10),count(10)
          real    f(n1)
          character*(*) fname
        end subroutine
      end interface
c
c -------------------------------------------------------

      print*,'HERE: read_unidata_ruc_hyb'

      istatus = 1

      call s_len(cdfname,slen)

      print*,'cdfname: ',cdfname(1:slen)

      call get_nvaltimes_unidata(cdfname,ntbg,ivaltimes,istatus)
c
      print*,'opening cdf file: ',cdfname(1:slen)

      rcode = NF_OPEN(cdfname,NF_NOWRITE,ncid)
      if(rcode.ne.NF_NOERR) then
         print *, NF_STRERROR(rcode)
         print *,'NF_OPEN ',cdfname(1:slen)
         return
      endif

      read(af,'(i4)') nn

      n=1
      do while(n.lt.ntbg.and.ivaltimes(n)/3600.ne. nn)
         n=n+1
      enddo
      if(ivaltimes(n)/3600.ne.nn) then

         print*,'ERROR: No record valid at requested time '
         print*,'ntbg/nn/af/n/ivaltimes(n) ',ntbg,' ',nn,' ',af,
     &' ',n,' ',ivaltimes(n)

         rcode= NF_CLOSE(ncid)
         if(rcode.ne.NF_NOERR) then
            print *, NF_STRERROR(rcode)
            print *,'in NF_GET_VAR: ',cmodel
            return
         endif

         goto 999

      else

         print*,'Found valid record at ivaltime'
         print*,'ntbg/nn/af/n/ivaltimes(n) ',ntbg,' ',nn,' ',af,
     &' ',n,' ',ivaltimes(n) 
         print*
      endif

      ! Get the pressure levels for this data
      nf_status = NF_INQ_VARID(ncid,'P_hybr',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'in level '
        istatus = 0
        return
      endif
      nf_status = NF_GET_VAR_REAL(ncid,nf_vid,prbght)
 
      ! Get index for 2m and 10m winds
      ind2m = 1
      ind10m = 2
  
      ! Set some indices
      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=1
      count(3)=nzbght
      start(4)=n
      count(4)=1

      print*,'read ht'
      cvar='Z_hybr'
      call read_netcdf_real(ncid,cvar,nxbg*nybg*count(3),ht,start
     +     ,count,rcode)
      if(rcode.ne.NF_NOERR) then
         if(rcode.gt.-61)then
            print *, NF_STRERROR(rcode)
            print *,'in NF_GET_VAR (Z): ',cmodel
         else
            print *,'Missing HT data detected: return'
         endif
         print*
         return
      endif

c
c ****** Statements to fill tp.
c
      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=1
      count(3)=nzbgtp
      start(4)=n
      count(4)=1

      print*,'read vptmp'
      cvar='vptmp_hybr'
      call read_netcdf_real(ncid,cvar,nxbg*nybg*count(3),tp,start
     +     ,count,rcode)

      if(rcode.ne.NF_NOERR) then
         if(rcode.gt.-61)then
            print *, NF_STRERROR(rcode)
            print *,'in NF_GET_VAR (t): ',cmodel
         else
            print *,'Missing T data detected: return'
         endif
         print*
         return
      endif

c
c ****** Statements to fill sh.                           
c
      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=1
      count(3)=nzbgsh
      start(4)=n
      count(4)=1

      print*,'read qv'
      cvar='hum_mix_hybr'
      call read_netcdf_real(ncid,cvar,nxbg*nybg*count(3),sh,start
     +     ,count,rcode)

      if(rcode.ne.NF_NOERR) then
         if(rcode.gt.-61)then
            print *, NF_STRERROR(rcode)
            print *,'in NF_GET_VAR (qv): ',cmodel
         else
            print *,'Missing qv data detected: return'
         endif
         print*
         return
      endif
c
c ****** Statements to fill uw. 
c
      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=1
      count(3)=nzbguv
      start(4)=n
      count(4)=1

      print*,'read uw'
      cvar='u_hybr'
      call read_netcdf_real(ncid,cvar,nxbg*nybg*count(3),uw,start
     +     ,count,rcode)

      if(rcode.ne.NF_NOERR) then
         if(rcode.gt.-61)then
            print *, NF_STRERROR(rcode)
            print *,'in NF_GET_VAR (uw): ',cmodel
         else
            print *,'Missing U data detected: return'
         endif
         print*
         return
      endif

c
c ****** Statements to fill vw.                           
c
      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=1
      count(3)=nzbguv
      start(4)=n
      count(4)=1

      print*,'read vw'
      cvar='v_hybr'
      call read_netcdf_real(ncid,cvar,nxbg*nybg*count(3),vw,start
     +     ,count,rcode)

      if(rcode.ne.NF_NOERR) then
         if(rcode.gt.-61)then
            print *, NF_STRERROR(rcode)
            print *,'in NF_GET_VAR (vw): ',cmodel
         else
            print *,'Missing V data detected: return'
         endif
         print*
         return
      endif
c
c ****** Statements to fill ww.
c
      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=1
      count(3)=nzbgww
      start(4)=n
      count(4)=1

      print*,'read omega'
      cvar='omega_hybr'
      call read_netcdf_real(ncid,cvar,nxbg*nybg*count(3),ww  
     +  ,start,count,rcode)

      if(rcode.ne.NF_NOERR) then
         if(rcode.gt.-61)then
            print *, NF_STRERROR(rcode)
            print *,'in NF_GET_VAR (ww): ',cmodel
         else
            print *,'Missing ww data detected: continue without'
            print *,'Filling ww with 0.0'
            ww(:,:,:) = 0.0
         endif
         print*
      endif

C   Get 2m T and RH, 10m U and V
      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=ind2m
      count(3)=1
      start(4)=n
      count(4)=1
      print*,'read t_2m'
      lcmpsfcq=.true.
      cvar='T_fhg'
      call read_netcdf_real(ncid,cvar,nxbg*nybg,tp_sfc,start
     +     ,count,rcode)
      if(rcode.ne.NF_NOERR) then
         if(rcode.gt.-61)then
            print *, NF_STRERROR(rcode)
            print *,'in NF_GET_VAR (tp_sfc): ',cmodel
         endif
      endif

c  No surface height...so use lowest hybrid level heigh
      ht_sfc = ht(:,:,1)

      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=ind2m
      count(3)=1
      start(4)=n
      count(4)=1
      print*,'read sh_sfc'
      lcmpsfcq=.true.
      cvar='hum_mix_fhg'
      call read_netcdf_real(ncid,cvar,nxbg*nybg,sh_sfc,start
     +     ,count,rcode)
      if(rcode.ne.NF_NOERR) then
         if(rcode.gt.-61)then
            print *, NF_STRERROR(rcode)
            print *,'in NF_GET_VAR (sh_sfc): ',cmodel
         endif
      endif

      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=ind10m
      count(3)=1
      start(4)=n
      count(4)=1
      print*,'read uw_sfc'
      lcmpsfcq=.true.
      cvar='u_fhg'
      call read_netcdf_real(ncid,cvar,nxbg*nybg,uw_sfc,start
     +     ,count,rcode)
      if(rcode.ne.NF_NOERR) then
         if(rcode.gt.-61)then
            print *, NF_STRERROR(rcode)
            print *,'in NF_GET_VAR (uw_sfc): ',cmodel
         endif
      endif

      start(1)=1
      count(1)=nxbg
      start(2)=1
      count(2)=nybg
      start(3)=ind10m
      count(3)=1
      start(4)=n
      count(4)=1
      print*,'read vw_sfc'
      lcmpsfcq=.true.
      cvar='v_fhg'
      call read_netcdf_real(ncid,cvar,nxbg*nybg,vw_sfc,start
     +     ,count,rcode)
      if(rcode.ne.NF_NOERR) then
         if(rcode.gt.-61)then
            print *, NF_STRERROR(rcode)
            print *,'in NF_GET_VAR (vw_sfc): ',cmodel
         endif
      endif

c
c get sfc pressure field (use lowest hybrid level)
c
c
      pr_sfc = prbght(:,:,1)

c get mslp (this field name differs from one model to the other)
c
      if(cmodel(1:3).eq.'RUC')THEN

         cvar='Pm_msl'
         call read_netcdf_real(ncid,cvar,nxbg*nybg,mslp
     +           ,start,count,rcode)
         if(rcode.ne.NF_NOERR) then
            if(rcode.gt.-61)then
               print *, NF_STRERROR(rcode)
               print *,'in NF_GET_VAR (emsp): ',cmodel
            else
               print*,'Error status returned from read_netcdf_real'
            endif
            print *,'Missing emsp data detected'
            print*
            return
         endif

      endif

c
c *** Close netcdf file.
c
      rcode= NF_CLOSE(ncid)
      if(rcode.ne.NF_NOERR) then
         print *, NF_STRERROR(rcode)
         print *,'in NF_GET_VAR: ',cmodel
         return
      endif
c
ccc      endif

c
c *** Fill ouput arrays.
c ***   convert theta-v to temperature
c ***   convert mixing ratio to specific humidity
c ***   convert 3d pressure level arrays to from Pa to mb
c ***   load remaining pressure arrays
c
c for laps-lgb

      call s_len(ctype,lent)

      print*,'ctype ',ctype(1:lent)
     
      do j=1,nybg
        do i=1,nxbg
           do k = 1, nzbght

c             Conver pressure to mb and fill additional arrays
              prbght(i,j,k) = prbght(i,j,k) * 0.01
              prbgsh(i,j,k) = prbght(i,j,k)
              prbguv(i,j,k) = prbght(i,j,k)
              prbgww(i,j,k) = prbght(i,j,k)

c             Convert theta-v into virtual temp
              tv = tp(i,j,k) * (prbgsh(i,j,k)*0.001)**rcp
C             Convert virtual temp into temp
              tp(i,j,k) = tv/(1.0+0.61*sh(i,j,k))

c             Convert mixrat into spechum
              sh(i,j,k) = sh(i,j,k)/(1.+sh(i,j,k))

           enddo
c          Also convert surface mixrat into spechum
           sh_sfc(i,j)=sh_sfc(i,j)/(1.+sh_sfc(i,j))
         enddo
       enddo

      ic = nxbg/2
      jc = nybg/2
      do k = 1, nzbght
        print 850, k,prbght(ic,jc,k),ht(ic,jc,k),tp(ic,jc,k), 
     +              sh(ic,jc,k),uw(ic,jc,k),vw(ic,jc,k),ww(ic,jc,k)
      enddo
! Honglig Jiang: change from F7.5 to F8.5. 11/27/2013
 850  FORMAT(I2,1x,F6.1,1x,F7.1,1x,F5.1,1x,F8.5,1x,F5.1,1x,F5.1,1x,
     +   F8.5)
    
      istatus =  1
      if(0.eq.1) then
 900     print*,'ERROR: bad dimension specified in netcdf file'
         print*, (count(i),i=1,4)
         istatus=-1
      endif
 999  return
      end
