      subroutine read_dgprep(bgmodel,cmodel,path,fname,af,nx,ny,nz
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
     .         ,icm
c
      integer  it
      integer  lun
      integer  iostat,iostatus
      integer  nclen

      logical  lopen,lext

      real*4 ht(nx,ny,nz)      !height (m)
     .      ,tp(nx,ny,nz)      !temperature (K)
     .      ,sh(nx,ny,nz)      !specific humidity (kg/kg) 
     .      ,uw(nx,ny,nz)      !u-wind (m/s)
     .      ,vw(nx,ny,nz)      !v-wind (m/s)
     .      ,pr(nx,ny,nz)      !pressures (mb)
     .      ,pw(nx,ny,nz)      !precip h2o for /public AVN
     .      ,prk(nz)

      real*4 ht_sfc(nx,ny)
     .      ,pr_sfc(nx,ny)
     .      ,td_sfc(nx,ny)     !for /public AVN this is RH @ 2m agl.
     .      ,tp_sfc(nx,ny)     !for /public AVN this is T @ 2m agl.
     .      ,uw_sfc(nx,ny)
     .      ,vw_sfc(nx,ny)
     .      ,mslp(nx,ny)
     .      ,accs(nx,ny)       !accum snow for /public AVN
     .      ,sfc_dummy(nx,ny)  !only used for /public AVN ... holds T @ sfc

      real*4 p_levels(nz,nvarsmax)

      double precision isoLevel(nz),reftime,valtime

      real*4 mrsat
      real*4 esat,xe
c     real*4 rp_init
      real*4 prsfc
      real*4 qsfc
      real*4 make_td
      real*4 make_ssh
      real*4 ssh2
      real*4 t_ref
      real*4 pcnt
      real*4 r_missing_data
      real*4 bogus_sh
c
      character*(*) path,cmodel
      character*9   fname
      character*4   af
      character*2   gproj
      character*13  fname13_to_FA_filename,
     .              cfname13,cFA_filename
      character*3   c3ext,  c3_FA_ext
      character*132 origin,model,nav,grid,version
      character*255 filename,fname_index
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
c Td or rh liq/ice phase temp thresh
c ---------------------------------
      t_ref=-47.0
      bogus_sh = 0.00001
c
      call get_r_missing_data(r_missing_data,istatus)

      call s_len(cmodel,nclen)

      lun=10
c
      if(bgmodel.eq.6.or.bgmodel.eq.8)then

         call s_len(path,l)
         filename=path(1:l)//'/'//fname//af
         call s_len(filename,l)

         if(cmodel(1:nclen).eq.'AVN_AFWA_DEGRIB')then
c
c *** Open and read data index file; and AFWA database thing.
c
           fname_index=filename(1:l)//'.index'

           call readindexfile(fname_index,nvarsmax,nz,nvars,nlevs
     +,p_levels,ivarcoord,ivarid,istatus)
           if(istatus.lt.1)goto 995

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
           print *,'Reading - ',filename(1:l)
           open(lun,file=filename(1:l),status='old',
     .          form='unformatted',err=990)
           rewind(lun)

           if(bgmodel.eq.6)then

              call read_avn(lun,nx,ny,nz,tp,uw,vw,ht,sh
     +,nvarsmax,nvars,nlevs,ivarcoord,ivarid
     +,ht_sfc,pr_sfc,td_sfc,tp_sfc,uw_sfc,vw_sfc,mslp
     +,istatus)

           elseif(bgmodel.eq.8)then

              call read_nogaps(lun,nx,ny,nz
     + ,nvarsmax,nvars,nlevs,ivarcoord,ivarid
     + ,ht,tp,sh,uw,vw,ht_sfc,pr_sfc,td_sfc,tp_sfc
     + ,uw_sfc,vw_sfc,mslp,istatus)


           endif

c        else
c
c eta ingest currently disabled. J. Smart (9-2-98)
c
c           call read_eta(lun,nx,ny,nz,tp,uw,vw,ht,sh
c    +,ht_sfc,pr_sfc,td_sfc,tp_sfc,uw_sfc,vw_sfc,mslp
c    +,istatus)

         elseif(cmodel(1:nclen).eq.'AVN_FSL_NETCDF')then

               call read_avn_netcdf(filename, nz, 1, nx, ny,
     +     version, ACCS, ht , ht_sfc, pr_sfc, mslp, pw, sh, tp,
     +     tp_sfc, sfc_dummy, uw, vw, td_sfc, isoLevel,
     +     reftime, valtime, grid, model, nav, origin, istatus)

               call qcmodel_sh(nx,ny,1,td_sfc)  !td_sfc actually = RH for AVN.

               call qcmodel_sh(nx,ny,nz,sh)     !sh actually = RH for AVN.

         endif

      elseif(bgmodel.eq.3)then

         cfname13=fname//af
         cFA_filename=fname13_to_FA_filename(cfname13)
         call s_len(path,l)
         filename=path(1:l)//'/'//cFA_filename
         call s_len(filename,l)

         inquire(file=filename,exist=lext,opened=lopen,number=lun)
         if(.not.lext)then
            print*,'File does not exist: ',filename(1:l)
            goto 990
         endif
         if(lopen)then
            print*,'File is already open: ',filename(1:l)
            goto 990
         endif

         print*,'open and read FA file: ',filename(1:l)
         open(lun,file=filename(1:l),status='old'
     +,IOSTAT=IOSTATUS,err=991)

         call read_fa(lun,filename                      ! I
     .               ,nx,ny,nz                          ! I
     .               ,r_missing_data                    ! I
     .               ,prk                               ! O
     .               ,ht,tp,sh,uw,vw                    ! O
     .               ,mslp                              ! O
     .               ,istatus)                          ! O

      endif
 
      if(istatus .ne. 1)then
         print*,'Error reading data: ',cmodel(1:nclen),
     +' ',filename(1:l)
         return
      endif
c
c *** Fill pressure array and convert rh to specific humidity. 
c *** Note: sh and td_sfc arrays contain rh from AVN (bgmodel=6)
c           or FA model (bgmodel=3).
c
      if(bgmodel.eq.6.or.bgmodel.eq.3)then

         if(cmodel(1:nclen).eq.'AVN_FSL_NETCDF')then
            do k=1,nz
              prk(k)=isoLevel(k)
            enddo
         endif

         print*,'convert rh to q - 3D'
         icm=0
         do k=1,nz
         do j=1,ny
         do i=1,nx

            pr(i,j,k)=prk(k)
            if(bgmodel.eq.3)pr(i,j,k)=pr(i,j,k)/100.
            if(sh(i,j,k).gt.0.0 .and. nint(sh(i,j,k)).le.100.)then
               if(sh(i,j,k).gt.100.)sh(i,j,k)=100.
               sh(i,j,k)=make_ssh(pr(i,j,k)
     .                        ,tp(i,j,k)-273.15
     .                        ,sh(i,j,k)/100.,t_ref)*0.001

            else
               icm=icm+1
               sh(i,j,k)=bogus_sh
            endif
               
c           it=tp(i,j,k)*100
c           it=min(45000,max(15000,it)) 
c           xe=esat(it)
c           mrsat=0.00622*xe/(prk(k)-xe)        !Assumes that rh units are %
c           sh(i,j,k)=sh(i,j,k)*mrsat           !rh --> mr
c           sh(i,j,k)=sh(i,j,k)/(1.+sh(i,j,k))  !mr --> sh

         enddo
         enddo
         enddo

         if(icm.gt.0)then
            pcnt=float(icm)/float(nx*ny*nz)
            print*,'WARNING: suspect 3d rh data (#/%): ',icm,pcnt
     &,' Bogus Q used for these points'

         if(bgmodel.eq.3)then
            do j=1,ny
            do i=1,nx
               mslp(i,j)=mslp(i,j)/100.   !hpa
            enddo
            enddo
         endif

         endif

         if(bgmodel.eq.6)then
 
            print*,'convert rh to Td - sfc: bgmodel: ',bgmodel
            icm=0
            do j=1,ny
            do i=1,nx

               if(td_sfc(i,j).gt.0.0 .and. td_sfc(i,j).le.100.)then
                  prsfc=pr_sfc(i,j)/100.
                  qsfc=make_ssh(prsfc,tp_sfc(i,j)-273.15,td_sfc(i,j)/100.
     &,t_ref)
                  td_sfc(i,j)=make_td(prsfc,tp_sfc(i,j)-273.15,qsfc
     &,t_ref)+273.15
               else
                  td_sfc(i,j)=make_td(pr_sfc(i,j)/100.,tp_sfc(i,j)-273.15
     &,bogus_sh,t_ref)+273.15
                  icm=icm+1
               endif

c           it=tp_sfc(i,j)*100
c           it=min(45000,max(15000,it))
c           xe=esat(it)
c           mrsat=0.00622*xe/(prsfc-xe)         !Assumes that rh units are %
c           td_sfc(i,j)=td_sfc(i,j)*mrsat             !rh --> mr
c           td_sfc(i,j)=td_sfc(i,j)/(1.+td_sfc(i,j))  !mr --> sh

            enddo
            enddo

            if(icm.gt.0)then
               pcnt=float(icm)/float(nx*ny)
               print*,'WARNING: suspect 2d rh data (#/%): ',icm,pcnt
     &,' Bogus Q used for these points'
            endif

         endif

      elseif(bgmodel.eq.8)then
c
c *** Convert Td to sh and fill 3D pressure array.
c
         icm=0
         print*,'Convert Td to q - 3D'
         do k=1,nz
         do j=1,ny
         do i=1,nx
            pr(i,j,k)=prk(k)
            if (sh(i,j,k) .gt. -99999.) then
               
               sh(i,j,k)=ssh2(prk(k),tp(i,j,k)-273.15,sh(i,j,k)-273.15
     &,t_ref)*0.001

c              it=sh(i,j,k)*100
c              it=min(45000,max(15000,it))
c              xe=esat(it)
c              sh(i,j,k)=0.622*xe/(pr(i,j,k)-xe)
c              sh(i,j,k)=sh(i,j,k)/(1.+sh(i,j,k))

            else
               sh(i,j,k)=bogus_sh
               icm=icm+1
            endif
         enddo
         enddo
         enddo
         if(icm.gt.0)then
            pcnt=float(icm)/float(nx*ny*nz)
            print*,'WARNING: missing 3d NOGAPS Td data ',icm,pcnt
         endif
c
c check for missing Td NOGAPS data.
c
         icm=0
         do j=1,ny
         do i=1,nx
            if(td_sfc(i,j).eq.-99999.)then
               icm=icm+1
               td_sfc(i,j)=r_missing_data
            endif
         enddo
         enddo
         if(icm.gt.0)then
            pcnt=float(icm)/float(nx*ny)
            print*,'WARNING: missing 2d NOGAPS Td data ',icm,pcnt
         endif


      endif
c
c *** Fill the common block variables.
c
      if (bgmodel .eq. 3)then
         gproj='LC'
         nx_lc=nx
         ny_lc=ny
         nz_lc=nz
         lat1=10.0
         lat2=40.0
         lon0_lc=+120.
         sw(1)=15.879
         sw(2)=+112.545
         ne(1)=32.384
         ne(2)=+131.172 
      elseif (bgmodel .eq. 6) then
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
      elseif (bgmodel.eq.8)then
         gproj='LL'
         nx_ll=nx
         ny_ll=ny
         nz_ll=nz
         lat0=-90.0
         lon0_ll=0.0
         dlat=1.0
         dlon=1.0
      endif
c
      istatus=1
      return
c
990   continue
      print *,'Error finding dgprep file.'
      return
995   print*,'Error reading model index file.',filename(1:l)
      return
991   IF (IOSTATUS .NE. 0)THEN
         PRINT *,'ERROR READING ',FILENAME(1:l),' IO status is', 
     &  IOSTATUS
      END IF

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

      istatus=0

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
c qc for model rh=0.0. 
c
      call qcmodel_sh(nx,ny,1,sh_sfc)

      call qcmodel_sh(nx,ny,nz,sh)

      istatus=1
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

      istatus=0

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

      istatus=1
      return

50    print*,'error during read'
      return
      end
C
C
C
      subroutine qcmodel_sh(nx,ny,nz,sh)

      implicit none

      integer i,j,k,nx,ny,nz

      real*4  sh(nx,ny,nz)

      do k=1,nz
      do j=1,ny
      do i=1,nx

         if(sh(i,j,k).le.0.0)then
            if( (i.gt.1.and.i.lt.nx) .and.
     .          (j.gt.1.and.j.lt.ny) )then

                 sh(i,j,k)=(sh(i+1,j,k)+sh(i-1,j,k)+
     .                      sh(i,j+1,k)+sh(i,j-1,k) )/4.0
            elseif(i.eq.1)then
                 sh(i,j,k)=sh(i+1,j,k)
            elseif(j.eq.1)then
                 sh(i,j,k)=sh(i,j+1,k)
            elseif(i.eq.nx)then
                 sh(i,j,k)=sh(i-1,j,k)
            elseif(j.eq.ny)then
                 sh(i,j,k)=sh(i,j-1,k)
            endif
         endif
      enddo
      enddo
      enddo

      return
      end
