      subroutine get_conus_dims(fname,nx,ny,nz)
c
      implicit none
c
      include 'netcdf.inc'
c
      integer nx,ny,nz, j, ncid
      integer vdims(10),start(10),count(10)
      integer ntp,nvdim,nvs,lenstr,ndsize,rcode
c
      character*(*) fname
      character*31  dummy
c_______________________________________________________________________________
c
c *** Open the netcdf file.
c
      rcode=NF_OPEN(fname,NF_NOWRITE,ncid)
c
c *** Statements to fill nx.
c
      call NCVINQ(ncid,16,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do j=1,nvdim
         call NCDINQ(ncid,vdims(j),dummy,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      enddo
      rcode=NF_GET_VARA_INT(ncid,16,start,count,nx)
c
c *** Statements to fill ny.
c
      call NCVINQ(ncid,17,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do j=1,nvdim
         call NCDINQ(ncid,vdims(j),dummy,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      enddo
      rcode=NF_GET_VARA_INT(ncid,17,start,count,ny)
c
c *** Statements to fill nz.
c
      call NCVINQ(ncid,18,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do j=1,nvdim
         call NCDINQ(ncid,vdims(j),dummy,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      enddo
      rcode=NF_GET_VARA_INT(ncid,18,start,count,nz)
c
      rcode=nf_close(ncid)
c

      return 
      end
c
c===============================================================================
c
      subroutine read_conus_nws(path,fname,af,nx,ny,nz,
     .                          pr,ht,tp,sh,uw,vw,gproj,istatus)
c
      implicit none
c
      include 'netcdf.inc'
c
      integer nx,ny,nz,len
      integer nvs, ndsize, ncid, ntp, nvdim, lenstr
c
      integer rcode
c
c *** Output arrays.
c
      real*4 pr(nx,ny,nz)
     .      ,ht(nx,ny,nz)
     .      ,tp(nx,ny,nz)
     .      ,sh(nx,ny,nz)
     .      ,uw(nx,ny,nz)
     .      ,vw(nx,ny,nz)
     .      ,prn(nz)
c
      real*4 lci(nx,ny),lcj(nx,ny),
     .       lat(nx,ny),lon(nx,ny),
     .       angle(nx,ny)
c
      integer start(10),count(10)
      integer vdims(10) 
      character*31 dummy
c
      integer i,j,k,istatus
c
      character*(*) path
      character*256 bgname
      character*9   fname
      character*4   af
      character*2   gproj
c
      real*4 msgflg
c
c *** Common block variables for Lambert-conformal grid.
c
      integer nx_lc,ny_lc,nz_lc  !No. of LC domain grid points
      real*4 lat1,lat2,lon0,     !Lambert-conformal std lat1, lat, lon
     .       sw(2),ne(2)         !SW lat, lon, NE lat, lon
      common /lcgrid/nx_lc,ny_lc,nz_lc,lat1,lat2,lon0,sw,ne
c_______________________________________________________________________________
c
      msgflg=1.e30
c
c *** Open the netcdf file.
c
      call s_len(path,len)
      bgname=path(1:len)//'/'//fname//af
      rcode=NF_OPEN(bgname,NF_NOWRITE,ncid)
c
c *** Read netcdf data.
c *** Statements to fill prn.
c
      call NCVINQ(ncid,1,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do j=1,nvdim
         call NCDINQ(ncid,vdims(j),dummy,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      enddo
      rcode=NF_GET_VARA_REAL(ncid,1,start,count,prn)
c
c *** Statements to fill ht.
c
      call NCVINQ(ncid,2,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do j=1,nvdim
         call NCDINQ(ncid,vdims(j),dummy,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      enddo
      rcode=NF_GET_VARA_REAL(ncid,2,start,count,ht)
c
c *** Statements to fill tp.
c
      call NCVINQ(ncid,3,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do j=1,nvdim
         call NCDINQ(ncid,vdims(j),dummy,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      enddo
      rcode=NF_GET_VARA_REAL(ncid,3,start,count,tp)
c
c *** Statements to fill sh.
c
      call NCVINQ(ncid,4,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do j=1,nvdim
         call NCDINQ(ncid,vdims(j),dummy,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      enddo
      rcode=NF_GET_VARA_REAL(ncid,4,start,count,sh)
c
c *** Statements to fill uw.
c
      call NCVINQ(ncid,5,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do j=1,nvdim
         call NCDINQ(ncid,vdims(j),dummy,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      enddo
      rcode=NF_GET_VARA_REAL(ncid,5,start,count,uw)
c
c *** Statements to fill vw.
c
      call NCVINQ(ncid,6,dummy,ntp,nvdim,vdims,nvs,rcode)
      lenstr=1
      do j=1,nvdim
         call NCDINQ(ncid,vdims(j),dummy,ndsize,rcode)
         lenstr=lenstr*ndsize
         start(j)=1
         count(j)=ndsize
      enddo
      rcode=NF_GET_VARA_REAL(ncid,6,start,count,vw)
c
c *** Close netcdf file.
c
      rcode=nf_close(ncid)
c
c *** Fill missing value flag (in netcdf file, missing = -99999.)
c
      do k=1,nz
      do j=1,ny
      do i=1,nx
         pr(i,j,k)=prn(k)
         if (ht(i,j,k) .lt. -10000.) ht(i,j,k)=msgflg
         if (tp(i,j,k) .lt. -10000.) tp(i,j,k)=msgflg
         if (sh(i,j,k) .lt. -10000.) sh(i,j,k)=msgflg
         if (uw(i,j,k) .lt. -10000.) uw(i,j,k)=msgflg
         if (vw(i,j,k) .lt. -10000.) vw(i,j,k)=msgflg
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
      istatus=1
      return
      end
