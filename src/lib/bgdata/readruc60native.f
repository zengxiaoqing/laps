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
      subroutine read_ruc60_native(path,fname,af,nx,ny,nz,
     .                   pr,ht,tp,sh,uw,vw,gproj,lon0_ps,istatus)
c
c *** Subroutine to read 60 km ruc data on the native polar stereographic,
c        hybrid-b grid.
c *** Code modified from B. Schwartz auto netcdf generator.
c
      implicit none
c

      include 'netcdf.inc'
      integer ncid, ntp, nvdim, lenstr, nvs, ndsize
c
      integer nx,ny,nz,rcode
      real   cp,g,cpog
      parameter (cp=1004.686,g=9.80665,cpog=cp/g)
c
c *** RUC arrays.
c
      real   pr(nx,ny,nz),       !Output ruc pressure (mb)
     .       ht(nx,ny,nz),       !Output ruc height (m)
     .       tp(nx,ny,nz),       !Output ruc temperature (K)
     .       sh(nx,ny,nz),       !Output ruc specific humidity (kg/kg)
     .       uw(nx,ny,nz),       !Output ruc u-wind (m/s)
     .       vw(nx,ny,nz),       !Output ruc v-wind (m/s)
     .       th(nx,ny,nz),       !ruc virtual potential temperature
     .       pc(nx,ny,nz)        !ruc condensation pressure
c
      real   psi(nx,ny),psj(nx,ny),
     .       lat(nx,ny),lon(nx,ny),
     .       mr,tv
c
      integer start(10),count(10)
      integer vdims(10)
      integer ndims,nvars,ngatts,recdim,nrecs
      character*31 dummy
c
      integer i,j,k,l,istatus
c
      character*(*) path
      character*9   fname
      character*4   af
      character*255 cdfname
      character*2   gproj
      
c
c *** Common block variables for polar-stereographic grid.
c
      integer nx_ps,ny_ps,nz_ps  !No. of PS domain grid points
      real   lat0,lon0,rota,       !Pol ste. std lat, lon and rotation
     .       sw(2),ne(2)           !SW lat, lon, NE lat, lon
      common /psgrid/nx_ps,ny_ps,nz_ps,lat0,lon0,rota,sw,ne
      real   lon0_ps             !returned for wind rotations
c_______________________________________________________________________________
c      
c *** Open the netcdf file.
c
c      l=index(path,' ')-1
      call s_len(path,l)
      cdfname=path(1:l)//'/'//fname//af
      print *,'Reading - ',cdfname(1:l+14)
      rcode=NF_OPEN(cdfname,NF_NOWRITE,ncid)
      call NCINQ(ncid,ndims,nvars,ngatts,recdim,rcode)
      call NCDINQ(ncid,recdim,dummy,nrecs,rcode)
      if (nrecs .lt. 1) then
         print *,'Not enough records in netcdf file.'
         istatus=0
         return
      endif
c
c *** Read netcdf data.
c *** Statements to fill uw.
c
      call NCVINQ(ncid,1,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do j=1,nvdim
         call NCDINQ(ncid,vdims(j),dummy,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      enddo
      rcode=NF_GET_VARA_REAL(ncid,1,start,count,uw)
c
c *** Statements to fill vw.
c
      call NCVINQ(ncid,2,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do j=1,nvdim
         call NCDINQ(ncid,vdims(j),dummy,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      enddo
      rcode=NF_GET_VARA_REAL(ncid,2,start,count,vw)
c
c *** Statements to fill ht (Mont. stream fucntion / g).
c
      call NCVINQ(ncid,3,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do j=1,nvdim
         call NCDINQ(ncid,vdims(j),dummy,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      enddo
      rcode=NF_GET_VARA_REAL(ncid,3,start,count,ht)
c
c *** Statements to fill pc.
c
      call NCVINQ(ncid,4,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do j=1,nvdim
         call NCDINQ(ncid,vdims(j),dummy,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      enddo
      rcode=NF_GET_VARA_REAL(ncid,4,start,count,pc)
c
c *** Statements to fill pr.
c
      call NCVINQ(ncid,5,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do j=1,nvdim
         call NCDINQ(ncid,vdims(j),dummy,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      enddo
      rcode=NF_GET_VARA_REAL(ncid,5,start,count,pr)
c
c *** Statements to fill th.
c
      call NCVINQ(ncid,6,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do j=1,nvdim
         call NCDINQ(ncid,vdims(j),dummy,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      enddo
      rcode=NF_GET_VARA_REAL(ncid,6,start,count,th)
c
c *** Close netcdf file.
c
      rcode= NF_CLOSE(ncid)
c
c *** Convert Pascals to mb.
c *** Compute temp and sh from thetav and pc.
c *** Compute height from msf.
c
      do k=1,nz
      do j=1,ny
      do i=1,nx
         pr(i,j,k)=pr(i,j,k)*0.01
         call thvpc2tq(th(i,j,k),pc(i,j,k),pr(i,j,k),
     .                 tp(i,j,k),sh(i,j,k))
         mr=sh(i,j,k)/(1.-sh(i,j,k))
         tv=tp(i,j,k)*(1.+0.61*mr)
         ht(i,j,k)=ht(i,j,k)-cpog*tv
      enddo
      enddo
      enddo
c
c *** Fill the Polar-Stereographic common block variables.
c
      gproj='PS'
      nx_ps=nx
      ny_ps=ny
      nz_ps=nz
      rota=0.
      lat0=90.0
      lon0=-105.0
      lon0_ps=lon0
      sw(1)=22.83730698
      sw(2)=-120.4905014
      ne(1)=45.98867416
      ne(2)=-60.82944107
c
c *** Convert ruc winds from grid north to true north.
c commented JS 11-09-00.
c     do j=1,ny
c     do i=1,nx
c        psi(i,j)=float(i)
c        psj(i,j)=float(j)
c     enddo
c     enddo
c     call psij_2_latlon(nx*ny,psi,psj,lat,lon)
c
c     call uvgrid_to_uvtrue_a(uw,vw,lon,lon0,nx,ny,nz)
c
      istatus=1
      return
      end
