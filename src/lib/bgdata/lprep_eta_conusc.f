      subroutine lprep_eta_conusc(nx,ny,nz,pr,tp,rh
     .         ,tp_sfc,pr_sfc,rh_sfc
     .         ,gproj,lon0_lc,lat1_lc,lat2_lc,istatus)
c
      implicit none
c
c     integer nx,ny,nz,nvars,recdim,i,j,k,istatus
      integer nx,ny,nz,i,j,k,istatus
c
      real*4 tp(nx,ny,nz),
     .       rh(nx,ny,nz),   ! in as rh out as sh
     .       pr(nx,ny,nz),
     .      tp_sfc(nx,ny),
     .      rh_sfc(nx,ny),   ! in as rh out as sh
     .      pr_sfc(nx,ny)

c     real make_ssh
c
c
      character*2   gproj
      include 'bgdata.inc'
      real xe, mrsat
      integer it
c
c *** Common block variables for Lambert-conformal grid.
c
      integer nx_lc,ny_lc,nz_lc  !No. of LC domain grid points
      real*4 lat1,lat2,lon0,       !Lambert-conformal std lat1, lat, lon
     .       sw(2),ne(2)           !SW lat, lon, NE lat, lon
      common /lcgrid/nx_lc,ny_lc,nz_lc,lat1,lat2,lon0,sw,ne
      real*4 lon0_lc             !returned for wind rotations
      real*4 lat1_lc,lat2_lc     !    "
c
c_______________________________________________________________________________
c  pr is input as a single verticle column output as 3d array
c
      istatus = 1

      do k=1,nz
         pr(1,1,k) = pr(k,1,1)
      enddo

c
c convert ua rh to sh
c
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
c               rh(i,j,k)=make_ssh(pr(i,j,k),tp(i,j,k)-273.15,
c     +              rh(i,j,k)/100.,0.0)*0.001
            enddo
         enddo
      enddo
c
c convert sfc rh to sh
c
      do j=1,ny
         do i=1,nx
               it=tp_sfc(i,j)*100
               it=min(45000,max(15000,it))
               xe=esat(it)
               mrsat=0.00622*xe/(pr_sfc(i,j)*0.01-xe)
               rh_sfc(i,j)=rh_sfc(i,j)*mrsat
               rh_sfc(i,j)=rh_sfc(i,j)/(1.+rh_sfc(i,j)) 
c               rh_sfc(i,j)=make_ssh(pr_sfc(i,j)/100.,tp_sfc(i,j)-273.15,
c     +              rh_sfc(i,j)/100.,0.0)*0.001
         enddo
      enddo
      

c
c *** Fill Lambert-conformal common block variables.
c
      gproj='LC'
      nx_lc=nx
      ny_lc=ny
      lat1=25.0
      lat1_lc=lat1
      lat2=25.0
      lat2_lc=lat2
      lon0=-95.0
      lon0_lc=lon0
      sw(1)=12.19
      sw(2)=-133.459
      ne(1)=57.29
      ne(2)=-49.3849
c
      istatus = 0
      return
      end




