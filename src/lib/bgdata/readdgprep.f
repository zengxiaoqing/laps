      subroutine read_dgprep(bgmodel,path,fname,af,nx,ny,nz
     .                      ,pr,ht,tp,sh,uw,vw
     .                      ,gproj,istatus)

c
      implicit none
c
      integer   bgmodel,nx,ny,nz
     .         ,i,j,k,l,istatus
c
      integer  it
      integer  lun

      real*4 ht(nx,ny,nz)      !height (m)
     .      ,tp(nx,ny,nz)      !temperature (K)
     .      ,sh(nx,ny,nz)      !specific humidity (kg/kg) 
     .      ,uw(nx,ny,nz)      !u-wind (m/s)
     .      ,vw(nx,ny,nz)      !v-wind (m/s)
     .      ,pr(nx,ny,nz)      !pressures (mb)
     .      ,prk(nz)

      real*4 dummy(nx,ny)

      real*4 mrsat
      real*4 esat,xe
      real*4 rp_init
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
      if(bgmodel.eq.6)then
         prk( 1)=1000.
         prk( 2)= 975.
         prk( 3)= 950.
         prk( 4)= 925.
         prk( 5)= 900.
         prk( 6)= 850.
         prk( 7)= 800.
         prk( 8)= 750.
         prk( 9)= 700.
         prk(10)= 650.
         prk(11)= 600.
         prk(12)= 550.
         prk(13)= 500.
         prk(14)= 450.
         prk(15)= 400.
         prk(16)= 350.
         prk(17)= 300.
         prk(18)= 250.
         prk(19)= 200.
         prk(20)= 150.
         prk(21)= 100.
         prk(22)=  70.
         prk(23)=  50.
         prk(24)=  30.
         prk(25)=  20.
         prk(26)=  10.
      else
         rp_init=1000.
         do k=1,nz
            prk(k)=rp_init-((k-1)*25.)
         enddo
      endif
c
c_______________________________________________________________________________
c
c *** Open data file.
c
      call s_len(path,l)
      filename=path(1:l)//'/'//fname//af
      l=l+14
      print *,'Reading - ',filename(1:l)
      lun=1
      open(lun,file=filename(1:l),status='old',
     .     form='unformatted',err=990)
      rewind(1)
c     read(1) nxf,nyf,nzf
c     if (nx .ne. nxf .or. ny .ne. nyf .or.
c    .    nz .ne. nzf ) then
c        print *,'Grid dimension mismatch.'
c        print *,'   LAPS BG  nx, ny, nz =',nx,ny,nz
c        print *,'   DGPREP   nx, ny, nz =',nxf,nyf,nzf
c        print *,'abort ...'
c        stop
c     endif
c     read(1) prk

      if(bgmodel.eq.6)then
         call read_avn(lun,nx,ny,nz,tp,uw,vw,ht,sh
     +,istatus)
      else
         call read_eta(lun,nx,ny,nz,tp,uw,vw,ht,sh
     +,istatus)
      endif
 
      if(istatus .eq. 0)then
         print*,'data not read properly'
         return
      endif
c
c *** Fill pressure array and
c *** Convert 3d rh to specific humidity.
c
      print*,'convert rh to sh'
      do k=1,nz-1
      do j=1,ny
      do i=1,nx
         pr(i,j,k)=prk(k)
         it=tp(i,j,k)*100
         it=min(45000,max(15000,it))
         xe=esat(it)
         mrsat=0.00622*xe/(prk(k)-xe)        !Assumes that rh units are %
         sh(i,j,k)=sh(i,j,k)*mrsat           !rh --> mr
         sh(i,j,k)=sh(i,j,k)/(1.+sh(i,j,k))  !mr --> sh
      enddo
      enddo
      enddo
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
      end
c
c ********************************************************
      subroutine read_avn(lun,nx,ny,nz,tp,uw,vw,ht,sh
     +,istatus)

      implicit none

      integer   nx,ny,nz
     .         ,i,j,k,istatus
c
      integer  lun
      real*4 ht(nx,ny,nz)      !height (m)
     .      ,tp(nx,ny,nz)      !temperature (K)
     .      ,sh(nx,ny,nz)      !specific humidity (kg/kg)
     .      ,uw(nx,ny,nz)      !u-wind (m/s)
     .      ,vw(nx,ny,nz)      !v-wind (m/s)

      real*4 dummy(nx,ny)

      istatus=1

      print*,'Read T'
      do k=nz,2,-1
         read(lun,err=50) ((tp(i,j,k-1),i=1,nx),j=1,ny)
      enddo
      read(lun,err=50) ((dummy(i,j),i=1,nx),j=1,ny)
      print*,'Read u'
      do k=nz,2,-1
         read(lun,err=50) ((uw(i,j,k-1),i=1,nx),j=1,ny)
      enddo
      read(lun,err=50) ((dummy(i,j),i=1,nx),j=1,ny)
      print*,'Read v'
      do k=nz,2,-1
         read(lun,err=50) ((vw(i,j,k-1),i=1,nx),j=1,ny)
      enddo
      read(lun,err=50) ((dummy(i,j),i=1,nx),j=1,ny)
      print*,'Read Height'
      do k=nz,2,-1
         read(lun,err=50) ((ht(i,j,k-1),i=1,nx),j=1,ny)
      enddo
      read(lun,err=50) ((dummy(i,j),i=1,nx),j=1,ny)
      print*,'Read RH'
      do k=nz-10,2,-1    !27-10=17 -> prk(17)=300mb = first moisture level.
         read(lun,err=50) ((sh(i,j,k-1),i=1,nx),j=1,ny)
      enddo
c
c As at AFWA, rh above level  (300mb)=10%
c
c     print*,'set upper level rh to 10%'
      do k=18,nz
      do j=1,ny
      do i=1,nx
         sh(i,j,k)=10.0
      enddo
      enddo
      enddo

      return

50    print*,'error during read'
      return
      end
c
c********************************************************
      subroutine read_eta(lun,nx,ny,nz,tp,uw,vw,ht,sh
     +,istatus)

      implicit none

      integer   nx,ny,nz
     .         ,i,j,k,istatus
c
      integer  lun
      real*4 ht(nx,ny,nz)      !height (m)
     .      ,tp(nx,ny,nz)      !temperature (K)
     .      ,sh(nx,ny,nz)      !specific humidity (kg/kg)
     .      ,uw(nx,ny,nz)      !u-wind (m/s)
     .      ,vw(nx,ny,nz)      !v-wind (m/s)

      istatus=1

      print*,'Read T'
      do k=nz,1,-1
         read(lun,err=50) ((tp(i,j,k),i=1,nx),j=1,ny)
      enddo
      print*,'Read u'
      do k=nz,1,-1
         read(lun,err=50) ((uw(i,j,k),i=1,nx),j=1,ny)
      enddo
      print*,'Read v'
      do k=nz,1,-1
         read(lun,err=50) ((vw(i,j,k),i=1,nx),j=1,ny)
      enddo
      print*,'Read Height'
      do k=nz,1,-1
         read(lun,err=50) ((ht(i,j,k),i=1,nx),j=1,ny)
      enddo
      print*,'Read RH'
      do k=nz,1,-1
         read(lun,err=50) ((sh(i,j,k),i=1,nx),j=1,ny)
      enddo

      return

50    print*,'error during read'
      istatus=0

      return
      end
