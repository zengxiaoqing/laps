      subroutine read_nogaps(bgmodel,path,fname,af,nx,ny,nz
     .                      ,pr,ht,tp,sh,uw,vw
     .                      ,ht_sfc,pr_sfc,sh_sfc,tp_sfc
     .                      ,uw_sfc,vw_sfc,mslp
     .                      ,gproj,istatus)

c
      implicit none
c
      integer nx,ny,nz,i,j,k,l,it,istatus
      integer bgmodel
c
c
      real*4 ht(nx,ny,nz),     !NOGAPS height (m)
     .       tp(nx,ny,nz),     !NOGAPS temperature (K)
     .       sh(nx,ny,nz),     !NOGAPS specific humidity (kg/kg) 
     .       uw(nx,ny,nz),     !NOGAPS u-wind (m/s)
     .       vw(nx,ny,nz),     !NOGAPS v-wind (m/s)
     .       pr(nx,ny,nz),     !NOGAPS pressures (mb)
     .       prk(nz)

      real*4 ht_sfc(nx,ny)
     .      ,pr_sfc(nx,ny)
     .      ,sh_sfc(nx,ny)
     .      ,tp_sfc(nx,ny)
     .      ,uw_sfc(nx,ny)
     .      ,vw_sfc(nx,ny)
     .      ,mslp(nx,ny)

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

c     the next line inserted by John Smart & Capt Bob Williams, 12/17/97
c     removed again J.Smart 6/14/98
c     read(16) (prk(k),k=1,nz)
c     end of insert

      print*,'Read 3-d variables'
      do k=1,nz
         read(16,err=50) ((tp(i,j,k),i=1,nx),j=1,ny)
      enddo

c     print*,'Read u'
      do k=1,nz
         read(16,err=50) ((uw(i,j,k),i=1,nx),j=1,ny)
      enddo

c     print*,'Read v'
      do k=1,nz
         read(16,err=50) ((vw(i,j,k),i=1,nx),j=1,ny)
      enddo

c     print*,'Read Td'
      do k=1,nz-9
         read(16,err=50) ((sh(i,j,k),i=1,nx),j=1,ny) !Read in as dew point.
      enddo
      do k=nz-8,nz
      do j=1,ny
      do i=1,nx
         sh(i,j,k)=-99999.
      enddo
      enddo
      enddo

c     print*,'Read ht'
      do k=1,nz
         read(16,err=50) ((ht(i,j,k),i=1,nx),j=1,ny)
      enddo
c
c read sfc avn variables
c
      print*,'read sfc variables'
      read(16,err=50) ((tp_sfc(i,j),i=1,nx),j=1,ny)
      read(16,err=50) ((uw_sfc(i,j),i=1,nx),j=1,ny)
      read(16,err=50) ((vw_sfc(i,j),i=1,nx),j=1,ny)
      read(16,err=50) ((ht_sfc(i,j),i=1,nx),j=1,ny)
      read(16,err=50) ((sh_sfc(i,j),i=1,nx),j=1,ny)
      read(16,err=50) ((mslp(i,j),i=1,nx),j=1,ny)

      close(16)
c
c *** Convert dew point to specific humidity.
c *** Fill pressure array.
c
      print*,'Convert Td to q'
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

      if(.false.)then
      do j=1,ny
      do i=1,nx
         pr_sfc(i,j)=pr_sfc(i,j)/100.
         it=sh_sfc(i,j)*100
         it=min(45000,max(15000,it))
         xe=esat(it)
         sh_sfc(i,j)=0.622*xe/(pr_sfc(i,j)-xe)
         sh_sfc(i,j)=sh_sfc(i,j)/(1.+sh_sfc(i,j))
      enddo
      enddo
      endif
c
c *** Fill the Lat-Lon common block variables.
c
      gproj='LL'
      nx_ll=nx
      ny_ll=ny
      nz_ll=nz
      lat0=-90.0
      lon0=0.0
      if(bgmodel.eq.3)then
         dlat=2.5
         dlon=2.5
      elseif(bgmodel.eq.8)then
         dlat=1.0
         dlon=1.0
      endif
c
      istatus=1
      print*
      print*

      return
c
50    print*,'Error reading nogaps files'
      istatus=0
      return

990   continue
      print *,'Error finding nogaps file.'
      istatus=0
c
      end
