cdis    Forecast Systems Laboratory
cdis    NOAA/OAR/ERL/FSL
cdis    325 Broadway
cdis    Boulder, CO     80303
cdis 
cdis    Forecast Research Division
cdis    Local Analysis and Prediction Branch
cdis    LAPS 
cdis 
cdis    This software and its documentation are in the public domain and 
cdis    are furnished "as is."  The United States government, its 
cdis    instrumentalities, officers, employees, and agents make no 
cdis    warranty, express or implied, as to the usefulness of the software 
cdis    and documentation for any purpose.  They assume no responsibility 
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis    
cdis    Permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  All modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making 
cdis    the modifications.  If significant modifications or enhancements 
cdis    are made to this software, the FSL Software Policy Manager  
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
      subroutine read_ruc40_native(path,fname,af,nx,ny,nz,
     .                             pr,ht,tp,sh,uw,vw,gproj,istatus)
c
c *** Subroutine to read 40 km ruc data on the native polar stereographic,
c        hybrid-b grid.
c *** Code modified from B. Schwartz auto netcdf generator.
c
      implicit none
c
      include 'netcdf.inc'
c
      integer*4 nx,ny,nz,rcode
      real*4 cp,g,r,cpog,kappa
      parameter (cp=1004.686,g=9.80665,r=287.053,cpog=cp/g,kappa=r/cp)
c
c *** RUC arrays.
c
      real*4 pr(nx,ny,nz),       !Output ruc pressure (mb)
     .       ht(nx,ny,nz),       !Output ruc height (m)
     .       tp(nx,ny,nz),       !Output ruc temperature (K)
     .       sh(nx,ny,nz),       !Output ruc specific humidity (kg/kg)
     .       uw(nx,ny,nz),       !Output ruc u-wind (m/s)
     .       vw(nx,ny,nz),       !Output ruc v-wind (m/s)
     .       th(nx,ny,nz)        !ruc virtual potential temperature
c
      real*4 psi(nx,ny),psj(nx,ny),
     .       lat(nx,ny),lon(nx,ny),
     .       angle(nx,ny),
     .       tv
c
      integer*4 start(10),count(10)
      integer vdims(10)
      integer ncid,ntp,nvdim,nvs,lenstr,ndsize
      integer ndims,nvars,ngatts,recdim,nrecs
      character*31 dummy
c
      integer*4 i,j,k,l,istatus, string_space
c
      character*(*) path
      character*9   fname
      character*4   af
      character*255 cdfname
      character*2   gproj
c
c *** Common block variables for Lambert-conformal grid.
c
      integer*4 nx_lc,ny_lc,nz_lc  !No. of LC domain grid points
      real*4 lat1,lat2,lon0,       !Lambert-conformal std lat1, lat, lon
     .       sw(2),ne(2)           !SW lat, lon, NE lat, lon
      common /lcgrid/nx_lc,ny_lc,nz_lc,lat1,lat2,lon0,sw,ne
c_______________________________________________________________________________
c      
c *** Open the netcdf file.
c
c      l=index(path,' ')-1
      l=string_space(path)-1
      cdfname=path(1:l)//'/'//fname//af
      print *,'Reading - ',cdfname(1:l+14)
      ncid=ncopn(cdfname,ncnowrit,rcode)
      call ncinq(ncid,ndims,nvars,ngatts,recdim,rcode)
      call ncdinq(ncid,recdim,dummy,nrecs,rcode)
      if (nrecs .lt. 1) then
         print *,'Not enough records in netcdf file.'
         istatus=0
         return
      endif
c
c *** Read netcdf data.
c *** Statements to fill uw.
c
      call ncvinq(ncid,1,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      enddo
      call ncvgt(ncid,1,start,count,uw,rcode)
c
c *** Statements to fill vw.
c
      call ncvinq(ncid,2,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      enddo
      call ncvgt(ncid,2,start,count,vw,rcode)
c
c *** Statements to fill ht (Mont. stream fucntion / g).
c
      call ncvinq(ncid,3,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      enddo
      call ncvgt(ncid,3,start,count,ht,rcode)
c
c *** Statements to fill pr.
c
      call ncvinq(ncid,4,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      enddo
      call ncvgt(ncid,4,start,count,pr,rcode)
c
c *** Statements to fill th.
c
      call ncvinq(ncid,5,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      enddo
      call ncvgt(ncid,5,start,count,th,rcode)
c
c *** Statements to fill sh (read in as mr).
c
      call ncvinq(ncid,6,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do j=1,nvdim
         call ncdinq(ncid,vdims(j),dummy,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      enddo
      call ncvgt(ncid,6,start,count,sh,rcode)
c
c *** Close netcdf file.
c
      call ncclos(ncid,rcode)
c
c *** Convert Pascals to mb.
c *** Compute tv from thetav.
c *** Compute height from msf.
c *** Compute tp from tv.
c *** Compute sh from mr.
c
      do k=1,nz
      do j=1,ny
      do i=1,nx
         pr(i,j,k)=pr(i,j,k)*0.01
         tv=th(i,j,k)*(pr(i,j,k)*0.001)**kappa
         ht(i,j,k)=ht(i,j,k)-cpog*tv
         tp(i,j,k)=tv/(1.+0.61*sh(i,j,k))
         sh(i,j,k)=sh(i,j,k)/(1.+sh(i,j,k))
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

      sw(1)=16.2810
      sw(2)=-126.1378
      ne(1)=55.4818
      ne(2)=-57.3794
c
c *** Convert ruc winds from grid north to true north.
c
      do j=1,ny
      do i=1,nx
         psi(i,j)=float(i)
         psj(i,j)=float(j)
      enddo
      enddo
      call psij_2_latlon(nx*ny,psi,psj,lat,lon)
c
      call uvgrid_to_uvtrue_a(uw,vw,lon,lon0,nx,ny,nz,angle)
c
      istatus=1
      return
      end
