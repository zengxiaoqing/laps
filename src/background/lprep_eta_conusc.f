      subroutine lprep_eta_conusc(nx,ny,nz,ht,pr,tp,uw,vw,rh,
     .                           gproj,istatus)
c
      implicit none
c
      integer nx,ny,nz,nvars,recdim,start(10),count(10),k
c
      real*4 ht(nx,ny,nz),
     .       tp(nx,ny,nz),
     .       uw(nx,ny,nz),
     .       vw(nx,ny,nz),
     .       rh(nx,ny,nz),   ! in as rh out as mr
     .       pr(nx,ny,nz),
     .      tmp(nx,ny,nz)
c
      integer vdims(10) !Allow up to 10 dimensions
      integer nvs,nvdim,ntp,ndsize,j,lenstr,ncid,
     .          nrecs,ngatts,ndims,ipr(nz), len, it, i, istatus
c
      character*2   gproj
      real*4 xe,mrsat
      include 'bgdata.inc'

c
c *** Common block variables for Lambert-conformal grid.
c
      integer nx_lc,ny_lc,nz_lc  !No. of LC domain grid points
      real*4 lat1,lat2,lon0,       !Lambert-conformal std lat1, lat, lon
     .       sw(2),ne(2)           !SW lat, lon, NE lat, lon
      common /lcgrid/nx_lc,ny_lc,nz_lc,lat1,lat2,lon0,sw,ne
c
c_______________________________________________________________________________
c  pr is input as a single verticle column output as 3d array
c
      do k=1,nz
         pr(1,1,k) = pr(k,1,1)
      enddo


      do k=1,nz
         do j=1,ny
            do i=1,nx
               pr(i,j,k)=pr(1,1,k)
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




