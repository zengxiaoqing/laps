      subroutine read_dgprep(bgmodel,path,fname,af,nx,ny,nz
     .                      ,pr,ht,tp,sh,uw,vw
     .                      ,gproj,istatus)

c
      implicit none
c
      integer   bgmodel,nx,ny,nz,nxf,nyf,nzf
     .         ,i,j,k,l,istatus
c
      integer  it
      real*4 ht(nx,ny,nz)      !height (m)
     .      ,tp(nx,ny,nz)      !temperature (K)
     .      ,sh(nx,ny,nz)      !specific humidity (kg/kg) 
     .      ,uw(nx,ny,nz)      !u-wind (m/s)
     .      ,vw(nx,ny,nz)      !v-wind (m/s)
     .      ,pr(nx,ny,nz)      !pressures (mb)
     .      ,prk(nz)

      real*4 xp
      real*4 mrsat
      real*4 esat,xe
c
      character*(*) path
      character*9   fname
      character*4   af
      character*2   gproj
      character*255 filename
c
c *** Common block variables for lat-lon grid.
c
      integer*4 nx_ll,ny_ll,nz_ll
      real*4 lat0,lon0_ll,dlat,dlon
      common /llgrid/nx_ll,ny_ll,nz_ll,lat0,lon0_ll,dlat,dlon
c
c *** Common block variables for lambert-conformal grid.
c
      integer*4 nx_lc,ny_lc,nz_lc
      real*4 lat1,lat2,lon0_lc,sw(2),ne(2)
      common /lcgrid/nx_lc,ny_lc,nz_lc,lat1,lat2,lon0_lc,sw,ne
c
      common /estab/esat(15000:45000)
c
      prk( 1)=1000.
      prk( 2)= 925.
      prk( 3)= 850.
      prk( 4)= 700.
      prk( 5)= 500.
      prk( 6)= 400.
      prk( 7)= 300.
      prk( 8)= 250.
      prk( 9)= 200.
      prk(10)= 150.
      prk(11)= 100.
      prk(12)=  70.
      prk(13)=  50.
      prk(14)=  30.
      prk(15)=  20.
      prk(16)=  10.

c_______________________________________________________________________________
c
c *** Open AVN file.
c
      call s_len(path,l)
      filename=path(1:l)//'/'//fname//af
      l=l+14
      print *,'Reading - ',filename(1:l)
      open(1,file=filename(1:l),status='old',
     .     form='unformatted',err=990)
      rewind(1)
c     read(1) nxf,nyf,nzf
      if (nx .ne. nxf .or. ny .ne. nyf .or.
     .    nz .ne. nzf) then
         print *,'Grid dimension mismatch.'
         print *,'   LAPS BG  nx, ny, nz =',nx,ny,nz
         print *,'   DGPREP   nx, ny, nz =',nxf,nyf,nzf
         print *,'Abort...'
c        stop
      endif
c     read(1) prk
      do k=nz,2,-1
         read(1) ((tp(i,j,k-1),i=1,nx),j=1,ny)
      enddo
      do k=nz,2,-1
         read(1) ((uw(i,j,k-1),i=1,nx),j=1,ny)
      enddo
      do k=nz,2,-1
         read(1) ((vw(i,j,k-1),i=1,nx),j=1,ny)
      enddo
      do k=nz,2,-1
         read(1) ((ht(i,j,k-1),i=1,nx),j=1,ny)
      enddo
      do k=nz,2,-1
         read(1) ((sh(i,j,k-1),i=1,nx),j=1,ny)
      enddo
      close(1)
c
c *** Fill pressure array and
c *** Convert 3d rh to specific humidity.
c
      if(bgmodel.eq.6)then

      do k=1,nz
      do j=1,ny
      do i=1,nx
         pr(i,j,k)=prk(k)
         it=tp(i,j,k)*100
         it=min(45000,max(15000,it))
         xe=esat(it)
         mrsat=0.00622*xe/(prk(k)-xe)         !Assumes that rh units are %
         sh(i,j,k)=sh(i,j,k)*mrsat           !rh --> mr
         sh(i,j,k)=sh(i,j,k)/(1.+sh(i,j,k))  !mr --> sh
      enddo
      enddo
      enddo

      endif
c
c *** Fill the common block variables.
c
      if (bgmodel .eq. 6) then
         gproj='LL'
         nx_ll=nx
         ny_ll=ny
         nz_ll=nz
         lat0=-90.0
         lon0_ll=0.0
         dlat=1.0
         dlon=1.0
      elseif (bgmodel .eq. 7) then
         gproj='LC'
         nx_lc=nx
         ny_lc=ny
         nz_lc=nz
         lat1=25.0
         lat2=25.0
         lon0_lc=-95.0
         sw(1)=12.19
         sw(2)=-133.459
         ne(1)=57.29
         ne(2)=-49.3849
      endif
c
      istatus=1
      return
c
990   continue
      print *,'Error finding dgprep file.'
      istatus=0
      return
c
      end
