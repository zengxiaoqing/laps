      subroutine get_sbn_model_id(filename,model,ivaltimes,ntbg)
      implicit none
      include 'netcdf.inc'
      character*(*) model, filename
      integer ivaltimes(*)
      integer nf_fid, nf_vid, nf_status, ntbg

C
C  Open netcdf File for reading
C
      nf_status = NF_OPEN(filename,NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'NF_OPEN rucsbn'
      endif

      nf_status=NF_INQ_UNLIMDIM(nf_fid,nf_vid)

      if(nf_status.ne.NF_NOERR) then
         print *, NF_STRERROR(nf_status)
      endif

      nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,ntbg)

      if(nf_status.ne.NF_NOERR) then
         print *, NF_STRERROR(nf_status)
      endif


      nf_status = NF_INQ_VARID(nf_fid,'model',nf_vid)
      if(nf_status.ne.NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *,'in var model'
      endif
      nf_status = NF_GET_VAR_TEXT(nf_fid,nf_vid,model)
      if(nf_status.ne.NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *,'in NF_GET_VAR_ model '
      endif
      nf_status=NF_INQ_VARID(nf_fid,'valtimeMINUSreftime',nf_vid)
      if(nf_status.ne.NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *,'in NF_GET_VAR_ model '
      endif
      nf_status=NF_GET_VARA_INT(nf_fid,nf_vid,1,ntbg,ivaltimes)
      if(nf_status.ne.NF_NOERR) then
         print *, NF_STRERROR(nf_status)
         print *,'in NF_GET_VAR_ model '
      endif
      
      

      nf_status = nf_close(nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'nf_close'
      endif

      return
      end

      subroutine get_sbn_dims
C
C  Fill all dimension values
C
C
C  Open netcdf File for reading
C
      nf_status = NF_OPEN(filename,NF_NOWRITE,nf_fid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'NF_OPEN rucsbn'
      endif

C
C Get size of levels_19
C
      nf_status = NF_INQ_DIMID(nf_fid,'levels_19',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim levels_19'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,levels_19)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim levels_19'
      endif
C
C Get size of levels_35
C
      nf_status = NF_INQ_DIMID(nf_fid,'levels_35',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim levels_35'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,levels_35)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim levels_35'
      endif
C
C Get size of levels_38
C
      nf_status = NF_INQ_DIMID(nf_fid,'levels_38',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim levels_38'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,levels_38)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim levels_38'
      endif
C
C Get size of levels_39
C
      nf_status = NF_INQ_DIMID(nf_fid,'levels_39',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim levels_39'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,levels_39)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim levels_39'
      endif
C
C Get size of levels_40
C
      nf_status = NF_INQ_DIMID(nf_fid,'levels_40',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim levels_40'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,levels_40)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim levels_40'
      endif
C
C Get size of record
C
      nf_status = NF_INQ_DIMID(nf_fid,'record',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim record'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,record)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim record'
      endif
C
C Get size of x
C
      nf_status = NF_INQ_DIMID(nf_fid,'x',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim x'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,x)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim x'
      endif
C
C Get size of y
C
      nf_status = NF_INQ_DIMID(nf_fid,'y',nf_vid)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim y'
      endif
      nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,y)
      if(nf_status.ne.NF_NOERR) then
        print *, NF_STRERROR(nf_status)
        print *,'dim y'
      endif
      return
      end
C
C

      subroutine read_conus_211(path,fname,af,nx,ny,nz,
     .                            nxbg,nybg,nzbg,ntbg,
     .                            pr,ht,tp,sh,uw,vw,gproj,istatus)

c
      implicit none
c
      include 'netcdf.inc'
      integer ncid, lenstr, ntp, nvdim, nvs, ndsize
c
      integer nx,ny,nz
      integer nxbg,nybg,nzbg(5),ntbg
c
      integer rcode
c
c *** Netcdf arrays.
c
cc      integer nxbg, nybg, nzbg1,nzbg2, nzbg3, nzbg4, ntbg
cc      parameter(nxbg=75,nybg=56,nzbg1=35,nzbg2=38,nzbg3=39)
cc      parameter(nzbg4=40,ntbg=5)

      real*4 htn(nxbg,nybg,nzbg(1)),
     .       rhn(nxbg,nybg,nzbg(2)),
     .       tpn(nxbg,nybg,nzbg(3)),
     .       uwn(nxbg,nybg,nzbg(4)),
     .       vwn(nxbg,nybg,nzbg(5)),
     .       prn(19)
c
      data prn/1000.,950.,900.,850.,800.,750.,700.,650.,600.,550.,
     .               500.,450.,400.,350.,300.,250.,200.,150.,100./
c
c *** Output arrays.
c
      real*4 pr(nx,ny,nz),
     .       ht(nx,ny,nz),
     .       tp(nx,ny,nz),
     .       sh(nx,ny,nz),
     .       uw(nx,ny,nz),
     .       vw(nx,ny,nz)
c
      real*4 lci(nx,ny),lcj(nx,ny),
     .       lat(nx,ny),lon(nx,ny),
     .       angle(nx,ny)
c
      integer start(10),count(10)
      integer vdims(10) 
      character*31 dummy
c
      integer i,j,k,n,ip,jp,ii,jj,kp1,it,istatus, slen
c
      character*(*) path
      character*4   af
      character*2   gproj
c
      real*4 msgflg
c
      real*4 xe,esat,mrsat
      common /estab/esat(15000:45000)
c
c *** Common block variables for Lambert-conformal grid.
c
      integer nx_lc,ny_lc,nz_lc  !No. of LC domain grid points
      real*4 lat1,lat2,lon0,       !Lambert-conformal std lat1, lat, lon
     .       sw(2),ne(2)           !SW lat, lon, NE lat, lon
      common /lcgrid/nx_lc,ny_lc,nz_lc,lat1,lat2,lon0,sw,ne
      character*13 fname13,fname9_to_wfo_fname13
      character*9 fname
      character*255 cdfname
c
ccc      save htn,tpn,rhn,uwn,vwn,prn,oldfname
c_______________________________________________________________________________
c
      msgflg=1.e30
c
c *** Open the netcdf file.
c
ccc      if (fname .ne. oldfname) then


C
C  Open netcdf File for reading
C
      fname13=fname9_to_wfo_fname13(fname)

      call s_len(path,slen)
      cdfname=path(1:slen)//'/'//fname13
      
      rcode = NF_OPEN(cdfname,NF_NOWRITE,ncid)
      if(rcode.ne.NF_NOERR) then
         print *, NF_STRERROR(rcode)
         print *,'NF_OPEN rucsbn'
      endif

c
c ****** Read netcdf data.
c ****** Statements to fill htn.
c
         read(af,'(i4)') n
         n=n/3+1
         if(n.gt.ntbg) return

         call NCVINQ(ncid,1,dummy,ntp,nvdim,vdims,nvs,rcode)
         lenstr=1
         do j=1,nvdim
            call NCDINQ(ncid,vdims(j),dummy,ndsize,rcode)
            lenstr=lenstr*ndsize
            start(j)=1
            count(j)=ndsize
         enddo
         
         if(count(1).ne.nxbg.or.count(2).ne.nybg.or.
     +      count(3).ne.nzbg(1).or.count(4).ne.ntbg) then
            goto 900
         endif         
         start(4)=n
         count(4)=1
      rcode=NF_GET_VARA_REAL(ncid,1,start,count,htn)
c
c ****** Statements to fill rhn.                           
c
         call NCVINQ(ncid,4,dummy,ntp,nvdim,vdims,nvs,rcode)
         lenstr=1
         do j=1,nvdim
            call NCDINQ(ncid,vdims(j),dummy,ndsize,rcode)
            lenstr=lenstr*ndsize
            start(j)=1
            count(j)=ndsize
         enddo
         if(count(1).ne.nxbg.or.count(2).ne.nybg.or.
     +      count(3).ne.nzbg(2).or.count(4).ne.ntbg) then
            goto 900
         endif         
         start(4)=n
         count(4)=1
      rcode=NF_GET_VARA_REAL(ncid,4,start,count,rhn)

c
c ****** Statements to fill tpn.                              
c
         call NCVINQ(ncid,7,dummy,ntp,nvdim,vdims,nvs,rcode)
         lenstr=1
         do j=1,nvdim
            call NCDINQ(ncid,vdims(j),dummy,ndsize,rcode)
            lenstr=lenstr*ndsize
            start(j)=1
            count(j)=ndsize
         enddo
         if(count(1).ne.nxbg.or.count(2).ne.nybg.or.
     +      count(3).ne.nzbg(3).or.count(4).ne.ntbg) then
            goto 900
         endif         
         start(4)=n
         count(4)=1
      rcode=NF_GET_VARA_REAL(ncid,7,start,count,tpn)
c
c ****** Statements to fill uwn.                           
c
         call NCVINQ(ncid,10,dummy,ntp,nvdim,vdims,nvs,rcode)
         lenstr=1
         do j=1,nvdim
            call NCDINQ(ncid,vdims(j),dummy,ndsize,rcode)
            lenstr=lenstr*ndsize
            start(j)=1
            count(j)=ndsize
         enddo
         if(count(1).ne.nxbg.or.count(2).ne.nybg.or.
     +      count(3).ne.nzbg(4).or.count(4).ne.ntbg) then
            goto 900
         endif         
         start(4)=n
         count(4)=1
      rcode=NF_GET_VARA_REAL(ncid,10,start,count,uwn)
c
c ****** Statements to fill vwn.                           
c
         call NCVINQ(ncid,13,dummy,ntp,nvdim,vdims,nvs,rcode)
         lenstr=1
         do j=1,nvdim
            call NCDINQ(ncid,vdims(j),dummy,ndsize,rcode)
            lenstr=lenstr*ndsize
            start(j)=1
            count(j)=ndsize
         enddo
         if(count(1).ne.nxbg.or.count(2).ne.nybg.or.
     +      count(3).ne.nzbg(5).or.count(4).ne.ntbg) then
            goto 900
         endif         
         start(4)=n
         count(4)=1
      rcode=NF_GET_VARA_REAL(ncid,13,start,count,vwn)
c
c *** Close netcdf file.
c
      rcode= NF_CLOSE(ncid)
c
ccc      endif

c
c *** Fill ouput arrays.
c *** Convert rh to sh.
c
      do k=1,nz
      do j=1,ny
      do i=1,nx
         pr(i,j,k)=prn(k)
         ht(i,j,k)=msgflg
         tp(i,j,k)=msgflg
         sh(i,j,k)=msgflg
         uw(i,j,k)=msgflg
         vw(i,j,k)=msgflg
      enddo
      enddo
      enddo
c         
c  For ruc the actual domain is smaller than the conus 211 projection
c  so nxbg and nybg are smaller than nx and ny
c
      if(nxbg.lt.nx.and.nybg.lt.ny) then
         ip=13
         jp=4 
      else
         ip=0
         jp=0
      endif

      n=1
      istatus=0
      do k=1,19
      do j=1,nybg
      do i=1,nxbg
         ii=i+ip
         jj=j+jp
         kp1=k+1
         if (htn(i,j,k) .gt. -1000. .and. 
     .       htn(i,j,k) .lt. 99999.) then
            ht(ii,jj,k)=htn(i,j,k)
            tp(ii,jj,k)=tpn(i,j,kp1)
            sh(ii,jj,k)=rhn(i,j,kp1)
            it=tp(ii,jj,k)*100
            it=min(45000,max(15000,it))
            xe=esat(it)
            mrsat=0.00622*xe/(prn(k)-xe)
            sh(ii,jj,k)=sh(ii,jj,k)*mrsat
            sh(ii,jj,k)=sh(ii,jj,k)/(1.+sh(ii,jj,k))
            uw(ii,jj,k)=uwn(i,j,kp1)
            vw(ii,jj,k)=vwn(i,j,kp1)
            istatus = 1
         endif
      enddo
      enddo
      enddo
      if(istatus .eq. 0) then
        print*, 'No valid data found for',fname, af
        return
      endif
cc      do jj=5,nybg
cc        do ii=14,nxbg
cc          if(uw(ii,jj,19).ge.msgflg .and. uw(ii,jj,18).lt.msgflg)
cc     +   then
cc            print*,'Filling top u level wind at ',ii,jj
cc            uw(ii,jj,19) = uw(ii,jj,18)
ccc          endif
cc          if(vw(ii,jj,19).ge.msgflg .and. vw(ii,jj,18).lt.msgflg)
cc     +   then
cc            print*,'Filling top v level wind at ',ii,jj
cc            vw(ii,jj,19) = vw(ii,jj,18)
cc          endif
cc
cc        enddo
cc      enddo
c
c *** Fill Lambert-conformal common block variables.
c
      gproj='LC'
      nx_lc=nx
      ny_lc=ny
      lat1=25.0
      lat2=25.0
      lon0=-95.0
      sw(1)=12.19
      sw(2)=-133.459
      ne(1)=57.29
      ne(2)=-49.3849
c
c *** Convert ruc winds from grid north to true north.
c
      istatus=1
      return

 900  print*,'ERROR: bad dimension specified in netcdf file'
      print*, (count(i),i=1,4)
      istatus=-1
      
      return
      end
