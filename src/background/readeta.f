      subroutine read_eta_conusc(bgpath,fname,af,nx,ny,nz,
     .                           pr, ht,tp,rh,uw,vw,gproj,istatus)
c
      implicit none
c
      integer*4 nx,ny,nz,nvars,recdim,start(10),count(10),k
c
      real*4 ht(nx,ny,nz),
     .       tp(nx,ny,nz),
     .       uw(nx,ny,nz),
     .       vw(nx,ny,nz),
     .       rh(nx,ny,nz),   ! in as rh out as mr
     .       pr(nx,ny,nz)
c
      integer vdims(10) !Allow up to 10 dimensions
      integer*4 nvs,nvdim,ntp,ndsize,j,lenstr,ncnowrit,ncopn,ncid,
     .          nrecs,ngatts,ndims,ipr(nz), len, it, i, istatus
c
      character*2   gproj
      character*132 etafile
      character*(*) bgpath
      character*9 fname
      character*4 af
      character*31  dummy
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
c_______________________________________________________________________________
c
c *** Open the netcdf file.
c
      call s_len(bgpath,len)
      etafile = bgpath(1:len)//'/'//fname//af

      ncid=ncopn(etafile,ncnowrit,istatus)
      call ncinq(ncid,ndims,nvars,ngatts,recdim,istatus)
      call ncdinq(ncid,recdim,dummy,nrecs,istatus)
      print *,'Reading - ',etafile
c
c *** Statements to fill uw.                            
c
      call ncvinq(ncid,10,dummy,ntp,nvdim,vdims,nvs,istatus)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),dummy,ndsize,istatus)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      enddo
      call ncvgt(ncid,10,start,count,uw,istatus)
c
c *** Statements to fill vw.                             
c
      call ncvinq(ncid,12,dummy,ntp,nvdim,vdims,nvs,istatus)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),dummy,ndsize,istatus)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      enddo
      call ncvgt(ncid,12,start,count,vw,istatus)
c
c *** Statements to fill ht.                             
c
      call ncvinq(ncid,2,dummy,ntp,nvdim,vdims,nvs,istatus)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),dummy,ndsize,istatus)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      enddo
      call ncvgt(ncid,2,start,count,ht,istatus)
c
c *** Statements to fill rh.                             
c
      call ncvinq(ncid,4,dummy,ntp,nvdim,vdims,nvs,istatus)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),dummy,ndsize,istatus)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      enddo
      call ncvgt(ncid,4,start,count,rh,istatus)
c
c *** Statements to fill tp.                             
c
      call ncvinq(ncid,7,dummy,ntp,nvdim,vdims,nvs,istatus)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),dummy,ndsize,istatus)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      enddo
      call ncvgt(ncid,7,start,count,tp,istatus)
c
c *** Statements to fill ipr.                       
c
      call ncvinq(ncid,31,dummy,ntp,nvdim,vdims,nvs,istatus)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),dummy,ndsize,istatus)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      enddo
      call ncvgt(ncid,31,start,count,ipr,istatus)
      
      do k=1,nz
         do j=1,ny
            do i=1,nx
               pr(i,j,k)=float(ipr(k))
               it=tp(i,j,k)*100
               it=min(45000,max(15000,it))
               xe=esat(it)
               mrsat=0.00622*xe/(pr(i,j,k)-xe)
               rh(i,j,k)=rh(i,j,k)*mrsat
               rh(i,j,k)=rh(i,j,k)/(1.+rh(i,j,k)) 
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
      istatus = 1
      return
      end
