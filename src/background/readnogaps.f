      subroutine read_nogaps(path,fname,af,nx,ny,nz,
     .                       pr,ht,tp,sh,uw,vw,
     .                       gproj,istatus)

c
      implicit none
c
      integer nx,ny,nz,i,j,k,l,it,istatus
c
c
      real*4 ht(nx,ny,nz),     !NOGAPS height (m)
     .       tp(nx,ny,nz),     !NOGAPS temperature (K)
     .       sh(nx,ny,nz),     !NOGAPS specific humidity (kg/kg) 
     .       uw(nx,ny,nz),     !NOGAPS u-wind (m/s)
     .       vw(nx,ny,nz),     !NOGAPS v-wind (m/s)
     .       pr(nx,ny,nz),     !NOGAPS pressures (mb)
     .       prk(nz)
c
      character*(*) path
      character*9   fname
      character*4   af
      character*2   gproj
      character*255 filename
c
      real*4 xe,esat
      common /estab/esat(15000:45000)
c
c *** Common block variables for lat-lon grid.
c
      integer nx_ll,ny_ll,nz_ll  !No. of LL domain grid points
      real*4 lat0,lon0,dlat,dlon   !Pol ste. std lat, lon and delta lat, lon
      common /llgrid/nx_ll,ny_ll,nz_ll,lat0,lon0,dlat,dlon
c_______________________________________________________________________________
c
c *** Fill NOGAPS pressure levels.
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
c
c *** Open nogaps file.
c
c      l=index(path//' ',' ')-1

      call s_len(path,l)
      filename=path(1:l)//'/'//fname//af
      l=l+14
      print *,'Reading - ',filename(1:l)
      open(16,file=filename(1:l),status='old',
     .     form='unformatted',err=990)
      rewind(16)
      read(16) (((ht(i,j,k),i=1,nx),j=1,ny),k=1,nz)
      read(16) (((tp(i,j,k),i=1,nx),j=1,ny),k=1,nz)
      read(16) (((sh(i,j,k),i=1,nx),j=1,ny),k=1,nz)  !Read in as dew point.
      read(16) (((uw(i,j,k),i=1,nx),j=1,ny),k=1,nz)
      read(16) (((vw(i,j,k),i=1,nx),j=1,ny),k=1,nz)
      close(1)
c
c *** Convert dew point to specific humidity.
c *** Fill pressure array.
c
      do k=1,nz
      do j=1,ny
      do i=1,nx
         pr(i,j,k)=prk(k)
         if (sh(i,j,k) .gt. -99999.) then
            it=sh(i,j,k)*100
            it=min(45000,max(15000,it))
            xe=esat(it)
            sh(i,j,k)=0.622*xe/(pr(i,j,k)-xe)
            sh(i,j,k)=sh(i,j,k)/(1.+sh(i,j,k))
         else
            if (pr(i,j,k) .lt. 300.) sh(i,j,k)=0.00001
         endif
      enddo
      enddo
      enddo
c
c *** Fill the Lat-Lon common block variables.
c
      gproj='LL'
      nx_ll=nx
      ny_ll=ny
      nz_ll=nz
      lat0=-90.0
      lon0=0.0
      dlat=2.5
      dlon=2.5
c
      istatus=1
      return
c
990   continue
      print *,'Error finding nogaps file.'
      istatus=0
c
      end
