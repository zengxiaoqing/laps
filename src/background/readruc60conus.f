      subroutine read_ruc60_conus(path,fname,af,nx,ny,nz,
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

c
      integer*4 rcode
c
c *** Netcdf arrays.
c
      real*4 htn(75,56,35,5),
     .       tpn(75,56,38,5),
     .       rhn(75,56,39,5),
     .       uwn(75,56,41,5),
     .       vwn(75,56,41,5),
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
      integer ndims,nvars,ngatts,recdim,nrecs
      character*31 dummy
c
      integer*4 i,j,k,l,n,ip,jp,ii,jj,kp1,it,istatus
c
      character*(*) path
      character*9   fname,oldfname
      character*4   af
      character*13  fname13,fname9_to_wfo_fname13
      character*255 cdfname
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
      save htn,tpn,rhn,uwn,vwn,prn,oldfname
c_______________________________________________________________________________
c
      msgflg=1.e30
c
c *** Open the netcdf file.
c
      if (fname .ne. oldfname) then
         fname13=fname9_to_wfo_fname13(fname)
         l=index(path,' ')-1
         cdfname=path(1:l)//'/'//fname13
         print *,'Reading - ',cdfname(1:l+14)
         ncid=ncopn(cdfname,ncnowrit,rcode)
         call ncinq(ncid,ndims,nvars,ngatts,recdim,rcode)
         call ncdinq(ncid,recdim,dummy,nrecs,rcode)
         if (nrecs .lt. 5) then
            print *,'Not enough records in netcdf file.'
            istatus=0
            return
         endif
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
         call ncvgt(ncid,13,start,count,vwn,rcode)
c
c *** Close netcdf file.
c
      call ncclos(ncid,rcode)
c
      endif
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
      ip=13
      jp=4 
      read(af,'(i4)') n
      n=n/3+1
      do k=1,19
      do j=1,56
      do i=1,75
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
         endif
      enddo
      enddo
      enddo
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
      do j=1,ny
      do i=1,nx
         lci(i,j)=float(i)
         lcj(i,j)=float(j)
      enddo
      enddo
      call lcij_2_latlon(nx*ny,lci,lcj,lat,lon)
c
      call uvgrid_to_uvtrue_a(uw,vw,lon,lon0,nx,ny,nz,angle)
c
      oldfname=fname
      istatus=1
      return
      end
