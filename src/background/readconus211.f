      subroutine open_sbn_netcdf(path, fname, ncid
     +                        , nrecs, istatus)  
      implicit none
      integer istatus
      include 'netcdf.inc'
      character*(*) path
      character*13 fname13,fname9_to_wfo_fname13
      character*9 fname
      character*255 cdfname
      character*31 dummy
      integer*4 ncid, slen, ndims, nvars, recdim, nrecs,ngatts
      integer*4 rcode

      if(len(path)+13.gt.len(cdfname)) then
        print*,'path too long in open_sbn_netcdf'
        stop
      endif

      fname13=fname9_to_wfo_fname13(fname)

      call s_len(path,slen)
      cdfname=path(1:slen)//'/'//fname13
      print *,'Reading - ',cdfname(1:slen+14)
      ncid=ncopn(cdfname,ncnowrit,rcode)
      call ncinq(ncid,ndims,nvars,ngatts,recdim,rcode)
      call ncdinq(ncid,recdim,dummy,nrecs,rcode)
      if (nrecs .le. 1) then
         print *,'Not enough records in netcdf file.'
         istatus=0
         call ncclos(ncid,istatus)
         return
      endif
      istatus = 1
      return
      end

      subroutine get_sbn_dims(ncid,nxbg,nybg,nzbg,ntbg,ivaltimes)
      integer*4 ncid,nxbg,nybg,nzbg(5),ntbg
      integer*4 id_fields(5), i, istat, vdims(10)
      integer ivaltimes(ntbg)      
      character*31 dummy
      data id_fields/1,4,7,10,13/
        
      i = ncvid(ncid,'valtimeMINUSreftime            ',istat)
      call ncvgt(ncid,i,1,ntbg,ivaltimes,istat)
      print *, ivaltimes
      do i=1,5
        call ncvinq(ncid,id_fields(i),dummy,ntp,nvdim,vdims,nvs,istat)

        call ncdinq(ncid,vdims(1),dummy,nxbg,rcode)
c        print *,'ndsize = ', ndsize
        call ncdinq(ncid,vdims(2),dummy,nybg,rcode)
c        print *,'ndsize = ', ndsize
        call ncdinq(ncid,vdims(3),dummy,nzbg(i),rcode)
c        print *,'ndsize = ', ndsize
cc        call ncdinq(ncid,vdims(4),dummy,ntbg,rcode)
c        print *,'ndsize = ', ndsize
        
c        print*, 'ntp = ',ntp
c        print*, 'nvdim = ',nvdim
c        print*, 'vdims = ',vdims
c        print*, 'nvs = ',nvs
      enddo
c      stop
      return 
      end


      subroutine read_conus_211(path,fname,af,nx,ny,nz,
     .                            nxbg,nybg,nzbg,ntbg,
     .                            pr,ht,tp,sh,uw,vw,gproj,istatus)

c
c *** Subroutine to read 60 km ruc data on the conus-c grid.
c *** Code modified from B. Schwartz auto netcdf generator.
c
      implicit none
c
      include 'netcdf.inc'
c
      integer*4 nx,ny,nz
      integer*4 nxbg,nybg,nzbg(5),ntbg
c
      integer*4 rcode
c
c *** Netcdf arrays.
c
cc      integer nxbg, nybg, nzbg1,nzbg2, nzbg3, nzbg4, ntbg
cc      parameter(nxbg=75,nybg=56,nzbg1=35,nzbg2=38,nzbg3=39)
cc      parameter(nzbg4=40,ntbg=5)

      real*4 htn(nxbg,nybg,nzbg(1),ntbg),
     .       tpn(nxbg,nybg,nzbg(2),ntbg),
     .       rhn(nxbg,nybg,nzbg(3),ntbg),
     .       uwn(nxbg,nybg,nzbg(4),ntbg),
     .       vwn(nxbg,nybg,nzbg(5),ntbg),
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
      integer*4 start(10),count(10)
      integer vdims(10) 
      integer ncid,ntp,nvdim,nvs,lenstr,ndsize
      character*31 dummy
c
      integer*4 i,j,k,n,ip,jp,ii,jj,kp1,it,istatus
c
      character*(*) path
      character*9   fname,oldfname
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
      integer*4 nx_lc,ny_lc,nz_lc  !No. of LC domain grid points
      real*4 lat1,lat2,lon0,       !Lambert-conformal std lat1, lat, lon
     .       sw(2),ne(2)           !SW lat, lon, NE lat, lon
      common /lcgrid/nx_lc,ny_lc,nz_lc,lat1,lat2,lon0,sw,ne
c
ccc      save htn,tpn,rhn,uwn,vwn,prn,oldfname
c_______________________________________________________________________________
c
      msgflg=1.e30
c
c *** Open the netcdf file.
c
ccc      if (fname .ne. oldfname) then

         call open_sbn_netcdf(path,fname,ncid,ntbg,rcode)
         if(rcode.ne.1) return
c
c ****** Read netcdf data.
c ****** Statements to fill htn.
c
         call ncvinq(ncid,1,dummy,ntp,nvdim,vdims,nvs,rcode)
         lenstr=1
         do j=1,nvdim
            call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
            lenstr=lenstr*ndsize
            start(j)=1
            count(j)=ndsize
         enddo
         if(count(1).ne.nxbg.or.count(2).ne.nybg.or.
     +      count(3).ne.nzbg(1).or.count(4).ne.ntbg) then
            goto 900
         endif         
         call ncvgt(ncid,1,start,count,htn,rcode)
c
c ****** Statements to fill rhn.                           
c
         call ncvinq(ncid,4,dummy,ntp,nvdim,vdims,nvs,rcode)
         lenstr=1
         do j=1,nvdim
            call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
            lenstr=lenstr*ndsize
            start(j)=1
            count(j)=ndsize
         enddo
         if(count(1).ne.nxbg.or.count(2).ne.nybg.or.
     +      count(3).ne.nzbg(2).or.count(4).ne.ntbg) then
            goto 900
         endif         
         call ncvgt(ncid,4,start,count,rhn,rcode)
c
c ****** Statements to fill tpn.                              
c
         call ncvinq(ncid,7,dummy,ntp,nvdim,vdims,nvs,rcode)
         lenstr=1
         do j=1,nvdim
            call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
            lenstr=lenstr*ndsize
            start(j)=1
            count(j)=ndsize
         enddo
         if(count(1).ne.nxbg.or.count(2).ne.nybg.or.
     +      count(3).ne.nzbg(3).or.count(4).ne.ntbg) then
            goto 900
         endif         
         call ncvgt(ncid,7,start,count,tpn,rcode)
c
c ****** Statements to fill uwn.                           
c
         call ncvinq(ncid,10,dummy,ntp,nvdim,vdims,nvs,rcode)
         lenstr=1
         do j=1,nvdim
            call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
            lenstr=lenstr*ndsize
            start(j)=1
            count(j)=ndsize
         enddo
         if(count(1).ne.nxbg.or.count(2).ne.nybg.or.
     +      count(3).ne.nzbg(4).or.count(4).ne.ntbg) then
            goto 900
         endif         

         call ncvgt(ncid,10,start,count,uwn,rcode)
c
c ****** Statements to fill vwn.                           
c
         call ncvinq(ncid,13,dummy,ntp,nvdim,vdims,nvs,rcode)
         lenstr=1
         do j=1,nvdim
            call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
            lenstr=lenstr*ndsize
            start(j)=1
            count(j)=ndsize
         enddo
         if(count(1).ne.nxbg.or.count(2).ne.nybg.or.
     +      count(3).ne.nzbg(5).or.count(4).ne.ntbg) then
            goto 900
         endif         
         call ncvgt(ncid,13,start,count,vwn,rcode)
c
c *** Close netcdf file.
c
      call ncclos(ncid,rcode)
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

      read(af,'(i4)') n
      n=n/3+1
      istatus=0
      do k=1,19
      do j=1,nybg
      do i=1,nxbg
         ii=i+ip
         jj=j+jp
         kp1=k+1
         if (htn(i,j,k,n) .gt. -1000. .and. 
     .       htn(i,j,k,n) .lt. 99999.) then
            ht(ii,jj,k)=htn(i,j,k,n)
            tp(ii,jj,k)=tpn(i,j,kp1,n)
            sh(ii,jj,k)=rhn(i,j,kp1,n)
            it=tp(ii,jj,k)*100
            it=min(45000,max(15000,it))
            xe=esat(it)
            mrsat=0.00622*xe/(prn(k)-xe)
            sh(ii,jj,k)=sh(ii,jj,k)*mrsat
            sh(ii,jj,k)=sh(ii,jj,k)/(1.+sh(ii,jj,k))
            uw(ii,jj,k)=uwn(i,j,kp1,n)
            vw(ii,jj,k)=vwn(i,j,kp1,n)
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
cc      do j=1,ny
cc      do i=1,nx
cc         lci(i,j)=float(i)
cc         lcj(i,j)=float(j)
cc      enddo
cc      enddo
cc      call lcij_2_latlon(nx*ny,lci,lcj,lat,lon)
c
cc      call uvgrid_to_uvtrue_a(uw,vw,lon,lon0,nx,ny,nz,angle)
c
cc      oldfname=fname
      istatus=1
      return
 900  print*,'ERROR: bad dimension specified in netcdf file'
      print*, (count(i),i=1,4)
      istatus=-1
      return
      end
