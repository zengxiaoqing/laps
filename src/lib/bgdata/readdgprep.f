      subroutine read_dgprep(bgmodel,path,fname,af,nx,ny,nz
     .                      ,pr,ht,tp,sh,uw,vw
     .                      ,ht_sfc,pr_sfc,td_sfc,tp_sfc
     .                      ,uw_sfc,vw_sfc,mslp
     .                      ,gproj,istatus)

c
      implicit none
c
      integer   nvarsmax
      parameter (nvarsmax=150)

      integer   ivarid(nvarsmax)
      integer   ivarcoord(nvarsmax)
      integer   nlevs(nvarsmax)
      integer   nvars

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

      real*4 ht_sfc(nx,ny)
     .      ,pr_sfc(nx,ny)
     .      ,td_sfc(nx,ny)
     .      ,tp_sfc(nx,ny)
     .      ,uw_sfc(nx,ny)
     .      ,vw_sfc(nx,ny)
     .      ,mslp(nx,ny)

      real*4 p_levels(nz,nvarsmax)

      real*4 mrsat
      real*4 esat,xe
c     real*4 rp_init
      real*4 prsfc
      real*4 qsfc
      real*4 make_td
      real*4 make_ssh
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
c reads model ".index" file. returns pressure of levels, variable id
c and number of levels for each model variable in file.  (J. Smart 7-6-98).
c
c_______________________________________________________________________________
c
c *** Open and read data index file.
c
      call s_len(path,l)
      filename=path(1:l)//'/'//fname//af//'.index'
      call s_len(filename,l)
      call readindexfile(filename,nvarsmax,nz,nvars,nlevs
     +,p_levels,ivarcoord,ivarid,istatus)
      if(istatus.ne.0)goto 995

      do j=1,nvars
         if(ivarid(j).eq.11.and.ivarcoord(j).eq.100)then
            do i=1,nlevs(j)
               prk(i)=p_levels(i,j)
            enddo
         endif
      enddo
c
c_______________________________________________________________________________
c
c *** Open and read data file.
c
      call s_len(path,l)
      filename=path(1:l)//'/'//fname//af
      l=l+14
      print *,'Reading - ',filename(1:l)
      lun=1
      open(lun,file=filename(1:l),status='old',
     .     form='unformatted',err=990)
      rewind(1)

      if(bgmodel.eq.6)then

         call read_avn(lun,nx,ny,nz,tp,uw,vw,ht,sh
     +,nvarsmax,nvars,nlevs,ivarcoord,ivarid
     +,ht_sfc,pr_sfc,td_sfc,tp_sfc,uw_sfc,vw_sfc,mslp
     +,istatus)

      elseif(bgmodel.eq.8.or.bgmodel.eq.3)then

         call read_nogaps(lun,nx,ny,nz
     + ,nvarsmax,nvars,nlevs,ivarcoord,ivarid
     + ,ht,tp,sh,uw,vw,ht_sfc,pr_sfc,td_sfc,tp_sfc
     + ,uw_sfc,vw_sfc,mslp,istatus)

c      else
c
c eta ingest currently disabled. J. Smart (9-2-98)
c
c        call read_eta(lun,nx,ny,nz,tp,uw,vw,ht,sh
c    +,ht_sfc,pr_sfc,td_sfc,tp_sfc,uw_sfc,vw_sfc,mslp
c    +,istatus)

      endif
 
      if(istatus .ne. 0)then
         print*,'data not read properly'
         return
      endif
c
c *** Fill pressure array and convert rh to specific humidity. 
c *** Note: sh and td_sfc arrays contain rh from AVN
c
      if(bgmodel.eq.6)then

         print*,'convert rh to q - 3D'
         do k=1,nz
         do j=1,ny
         do i=1,nx

            sh(i,j,k)=make_ssh(prk(k),
     .                         tp(i,j,k)-273.15,
     .                         sh(i,j,k)/100.,0.)*0.001

            pr(i,j,k)=prk(k)

c           it=tp(i,j,k)*100
c           it=min(45000,max(15000,it)) 
c           xe=esat(it)
c           mrsat=0.00622*xe/(prk(k)-xe)        !Assumes that rh units are %
c           sh(i,j,k)=sh(i,j,k)*mrsat           !rh --> mr
c           sh(i,j,k)=sh(i,j,k)/(1.+sh(i,j,k))  !mr --> sh

         enddo
         enddo
         enddo
 
         print*,'convert rh to Td - sfc'
         do j=1,ny
         do i=1,nx

            prsfc=pr_sfc(i,j)/100.
            qsfc=make_ssh(prsfc,tp_sfc(i,j)-273.15,td_sfc(i,j)/100.,0.)
            td_sfc(i,j)=make_td(prsfc,tp_sfc(i,j)-273.15,qsfc,0.)+273.15

c           it=tp_sfc(i,j)*100
c           it=min(45000,max(15000,it))
c           xe=esat(it)
c           mrsat=0.00622*xe/(prsfc-xe)         !Assumes that rh units are %
c           td_sfc(i,j)=td_sfc(i,j)*mrsat             !rh --> mr
c           td_sfc(i,j)=td_sfc(i,j)/(1.+td_sfc(i,j))  !mr --> sh

         enddo
         enddo

      elseif(bgmodel.eq.8.or.bgmodel.eq.3)then
c
c *** Convert Td to sh and fill 3D pressure array.
c
         print*,'Convert Td to q - 3D'
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
      elseif (bgmodel.eq.3.or.bgmodel.eq.8)then
         gproj='LL'
         nx_ll=nx
         ny_ll=ny
         nz_ll=nz
         lat0=-90.0
         lon0_ll=0.0
         if(bgmodel.eq.8)then
            dlat=1.0
            dlon=1.0
         else
            dlat=2.5
            dlon=2.5
         endif
      endif
c
      istatus=1
      return
c
990   continue
      print *,'Error finding dgprep file.'
      istatus=0
      return
995   print*,'Error reading model index file.',filename(1:l)
      istatus=0
      return
      end
c
c ********************************************************
      subroutine read_avn(lun,nx,ny,nz,tp,uw,vw,ht,sh
     +,nvarsmax,nvars,nlevs,ivarcoord,ivarid
     +,ht_sfc,pr_sfc,sh_sfc,tp_sfc,uw_sfc,vw_sfc,mslp
     +,istatus)

      implicit none

      integer   nx,ny,nz
     .         ,i,j,k,l,istatus
     .         ,nvarsmax,nvars
     .         ,nshl
     .         ,nlevs(nvarsmax)
     .         ,ivarcoord(nvarsmax)
     .         ,ivarid(nvarsmax)

c
      integer  lun

      real*4 ht(nx,ny,nz)      !height (m)
     .      ,tp(nx,ny,nz)      !temperature (K)
     .      ,sh(nx,ny,nz)      !specific humidity (kg/kg)
     .      ,uw(nx,ny,nz)      !u-wind (m/s)
     .      ,vw(nx,ny,nz)      !v-wind (m/s)

      real*4 ht_sfc(nx,ny)
     .      ,pr_sfc(nx,ny)
     .      ,sh_sfc(nx,ny)
     .      ,tp_sfc(nx,ny)
     .      ,uw_sfc(nx,ny)
     .      ,vw_sfc(nx,ny)
     .      ,mslp(nx,ny)

      real*4 dummy(nx,ny,nz)

      istatus=1

      print*,'read 3-d variables'
c nvar = 1
      do k=1,nz
         read(lun,err=50) ((tp(i,j,k),i=1,nx),j=ny,1,-1)
      enddo
c     read(lun,err=50) ((dummy(i,j),i=1,nx),j=1,ny)
c     print*,'Read u'
c = 2
      do k=1,nz
         read(lun,err=50) ((uw(i,j,k),i=1,nx),j=ny,1,-1)
      enddo
c     read(lun,err=50) ((dummy(i,j),i=1,nx),j=1,ny)
c     print*,'Read v'
c = 3
      do k=1,nz
         read(lun,err=50) ((vw(i,j,k),i=1,nx),j=ny,1,-1)
      enddo
c     read(lun,err=50) ((dummy(i,j),i=1,nx),j=1,ny)
c     print*,'Read Height'
c = 4
      do k=1,nz
         read(lun,err=50) ((ht(i,j,k),i=1,nx),j=ny,1,-1)
      enddo
c     read(lun,err=50) ((dummy(i,j),i=1,nx),j=1,ny)
c     print*,'Read RH'
c = 5
      nshl=nlevs(5)
      do k=1,nshl    ! -> prk(17)=300mb = last moisture level.
         read(lun,err=50) ((sh(i,j,k),i=1,nx),j=ny,1,-1)
      enddo
c
c read sfc avn variables
c
c = 6,7,8,9,10,11
      print*,'read sfc variables'
      read(lun,err=50) ((tp_sfc(i,j),i=1,nx),j=ny,1,-1)
      read(lun,err=50) ((uw_sfc(i,j),i=1,nx),j=ny,1,-1)
      read(lun,err=50) ((vw_sfc(i,j),i=1,nx),j=ny,1,-1)
      read(lun,err=50) ((ht_sfc(i,j),i=1,nx),j=ny,1,-1)
      read(lun,err=50) ((sh_sfc(i,j),i=1,nx),j=ny,1,-1)
      read(lun,err=50) ((mslp(i,j),i=1,nx),j=ny,1,-1)
c nvar = 12
c
      do l=12,nvars
        if(ivarid(l).eq.1.and.ivarcoord(l).eq.1)then
           read(lun,err=50) ((pr_sfc(i,j),i=1,nx),j=ny,1,-1)
           goto  188
        else
           do k=1,nlevs(l)
              read(lun,err=50)((dummy(i,j,k),i=1,nx),j=1,ny)
           enddo
        endif
      enddo
      print*,'Did not find mslp data!'

188   continue
c
c As at AFWA, rh above level  (300mb)=10%
c
c     print*,'set upper level rh to 10%'
      do k=nshl+1,nz
      do j=1,ny
      do i=1,nx
         sh(i,j,k)=10.0
      enddo
      enddo
      enddo

      istatus=0
      return

50    print*,'error during read'
      return
      end
c
c********************************************************
      subroutine read_eta(lun,nx,ny,nz,tp,uw,vw,ht,sh
     +,ht_sfc,pr_sfc,sh_sfc,tp_sfc,uw_sfc,vw_sfc,mslp
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

      real*4 ht_sfc(nx,ny)
     .      ,pr_sfc(nx,ny)
     .      ,sh_sfc(nx,ny)
     .      ,tp_sfc(nx,ny)
     .      ,uw_sfc(nx,ny)
     .      ,vw_sfc(nx,ny)
     .      ,mslp(nx,ny)

      istatus=1

      print*,'read 3-d variables'
      do k=1,nz
         read(lun,err=50) ((tp(i,j,k),i=1,nx),j=1,ny)
      enddo
c     print*,'Read u'
      do k=1,nz
         read(lun,err=50) ((uw(i,j,k),i=1,nx),j=1,ny)
      enddo
c     print*,'Read v'
      do k=1,nz
         read(lun,err=50) ((vw(i,j,k),i=1,nx),j=1,ny)
      enddo
c     print*,'Read Height'
      do k=1,nz
         read(lun,err=50) ((ht(i,j,k),i=1,nx),j=1,ny)
      enddo
c     print*,'Read RH'
      do k=1,nz
         read(lun,err=50) ((sh(i,j,k),i=1,nx),j=1,ny)
      enddo

c
c read eta sfc variables
c
      print*,'read sfc variables'
      read(lun,err=50) ((tp_sfc(i,j),i=1,nx),j=1,ny)
      read(lun,err=50) ((uw_sfc(i,j),i=1,nx),j=1,ny)
      read(lun,err=50) ((vw_sfc(i,j),i=1,nx),j=1,ny)
      read(lun,err=50) ((ht_sfc(i,j),i=1,nx),j=1,ny)
      read(lun,err=50) ((sh_sfc(i,j),i=1,nx),j=1,ny)
      read(lun,err=50) ((mslp(i,j),i=1,nx),j=1,ny)

      istatus=0
      return

50    print*,'error during read'
      return
      end
