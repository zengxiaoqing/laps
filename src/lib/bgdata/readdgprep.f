      subroutine read_dgprep(bgmodel,cmodel,path,fname,af,nx,ny,nz
     .                      ,pr,ht,tp,sh,uw,vw,ww
     .                      ,ht_sfc,pr_sfc,td_sfc,tp_sfc
     .                      ,uw_sfc,vw_sfc,mslp,istatus)

c
      implicit none
c
      integer   nvarsmax
      parameter (nvarsmax=150)

      integer   ivarid(nvarsmax)
      integer   ivarcoord(nvarsmax)
      integer   nlevs(nvarsmax)
      integer   nvars

      integer   bgmodel,nx,ny,nz,nzsh
     .         ,i,j,k,l,n,istatus,ksh
     .         ,icm(nz),icm_sfc
     .         ,icn3d
c
      real*4   shsum(nz),sumtot,sumsfctd
      real*4   shavg
      integer  it
      integer  lun
      integer  iostat,iostatus
      integer  nclen,nflen,lenc

      logical  lopen,lext

      real*4 ht(nx,ny,nz)      !height (m)
     .      ,tp(nx,ny,nz)      !temperature (K)
     .      ,sh(nx,ny,nz)      !specific humidity (kg/kg) 
     .      ,uw(nx,ny,nz)      !u-wind (m/s)
     .      ,vw(nx,ny,nz)      !v-wind (m/s)
     .      ,ww(nx,ny,nz)      !w-wind (pa/s)  !currently not available for NOGAPS
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
      real*4 r_bogus_sh
      real*4 ssh2
      real*4 t_ref
      real*4 pcnt
      real*4 r_missing_data
      real*4 rfill
      real*4 rmx2d,rmn2d
      integer imx,jmx,imn,jmn

c
      character*256 path
      character*132 cmodel
      character*200 fname
      character*10  cfname10,a9_to_yr_a10_time
      character*4   af
      character*2   gproj
 
      character*16  cFA_filename
      character*3   c3ext,  c3_FA_ext
      character*2   cwb_type
      character*132 origin,model,nav,grid,version
      character*255 filename,fname_index
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
      t_ref=-132.0
      rfill = -99999.
      r_bogus_sh = 1.0e-6
      istatus = 1
c
      call get_r_missing_data(r_missing_data,istatus)

      call s_len(cmodel,nclen)
      call s_len(fname,nflen)

      lun=10
c
      if(bgmodel.eq.6.or.bgmodel.eq.8)then

         call s_len(path,l)
         filename=path(1:l)//'/'//fname(1:nflen)   !//af
         call s_len(filename,l)

         if(cmodel(1:nclen).eq.'AVN_AFWA_DEGRIB'.or.
     +      cmodel(1:nclen).eq.'NOGAPS_AFWA_DEGRIB')then
c
c *** Open and read data index file; and AFWA database thing.
c
           fname_index=filename(1:l)//'.index'

           call readindexfile(fname_index,nvarsmax,nz,nvars,nlevs
     +,p_levels,ivarcoord,ivarid,istatus)
           if(istatus.eq.1)goto 995

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
           print *,'Opening ',filename(1:l)
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
     +     tp_sfc, sfc_dummy, uw, vw, ww, td_sfc, isoLevel,
     +     reftime, valtime, grid, model, nav, origin, istatus)

               nzsh=nz-5

               call qcmodel_sh(nx,ny,1,td_sfc)  !td_sfc actually = RH for AVN.

               call qcmodel_sh(nx,ny,nzsh,sh)     !sh actually = RH for AVN.

         endif

      elseif(bgmodel.eq.3)then

c note: library function fname13_to_FA_filename could be used
c       to covert the FA filename but currently is not.  J.Smart

         call downcase(cmodel(nclen-1:nclen),cwb_type)
         c3ext=c3_FA_ext(af) 
         cfname10=a9_to_yr_a10_time(fname,istatus)
         cFA_filename=cwb_type//cfname10(1:8)//fname(6:7)//'.'//c3ext
         call s_len(path,l)
         filename=path(1:l)//'/'//cFA_filename
         call s_len(filename,l)

         print*,'Opening FA file: ',filename(1:l)
         open(lun,file=filename(1:l),status='old'
     +,IOSTAT=IOSTATUS,err=990)

         if(cwb_type .eq. 're')then
            call read_fa(lun,filename                   ! I
     .               ,nx,ny,nz                          ! I
     .               ,r_missing_data                    ! I
     .               ,prk                               ! O
     .               ,ht,tp,sh,uw,vw                    ! O
     .               ,mslp                              ! O
     .               ,istatus)                          ! O

            call qcmodel_sh(nx,ny,nz,sh)
         elseif (cwb_type.eq.'nf')then
            print *,'read nf model with read_fa_nf' 
            call read_fa_nf(lun,filename,              ! I
     .                nx, ny, nz,                      ! I
     .                r_missing_data,                  ! I
     .                prk,                             ! O
     .                ht, tp, sh, uw, vw, ww,          ! O
     .                ht_sfc, tp_sfc, td_sfc,          ! O
     .                uw_sfc, vw_sfc, mslp,            ! O
     .                istatus )                        ! O
            call qcmodel_sh(nx,ny,nz,sh)
c           do k=1,nz
c              call get_mxmn_2d(nx,ny,ww,rmx2d,rmn2d
c    &,imx,jmx,imn,jmn)
c              print*,'k mx/mn ww ',k,rmx2d,rmn2d
c           enddo
         endif

         nzsh=nz

      endif

      if(istatus .eq. 1)then
         print*,'Error reading dgprep data file: ',cmodel(1:nclen),
     +' ',filename(1:l)
         return
      endif

c qc 

      do k=1,nz
         if(k.le.nzsh)then
            ksh=k
         else
            ksh=nzsh
         endif
         do j=1,ny
         do i=1,nx
            if((abs(ht(i,j,k))  .gt. 100000.) .or.
     +         (abs(tp(i,j,k))  .gt. 1000.)   .or.
     +             (tp(i,j,k)   .le. 0.)      .or.
     +         (abs(sh(i,j,ksh))  .ge. 101.)    .or.   !rh
     +              sh(i,j,ksh)   .eq. rfill    .or.
     +         (abs(uw(i,j,k))  .gt. 150.)    .or.
     +         (abs(vw(i,j,k))  .gt. 150.)        )then
 
               print*,'ERROR: Missing or bad value detected: ',i,j,k
               print*,'ht/tp/sh/uw/vw/ww: ',ht(i,j,k),tp(i,j,k)
     +                  ,sh(i,j,k), uw(i,j,k),vw(i,j,k),ww(i,j,k)

               istatus = 1
               return
            endif
         enddo
         enddo
      enddo

c qc - for rh.

c     do k=1,nzsh
c        do j=1,ny
c        do i=1,nx
c           if((abs(sh(i,j,k)) .ge. 101.)then
c              print*,'ERROR: Missing RH value detected: ',i,j,k
c              print*,'rh: ',sh(i,j,k)
c              istatus = 1
c              return
c           endif
c        enddo
c        enddo
c     enddo

 
      if(istatus .eq. 1)then
         print*,'Error reading dgprep data file: ',cmodel(1:nclen),
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
         do k=1,nz

          icm(k)=0
          shsum(k)=0.0
          do j=1,ny
          do i=1,nx

            pr(i,j,k)=prk(k)
            if(bgmodel.eq.3)pr(i,j,k)=pr(i,j,k)/100.

            if(sh(i,j,k) .gt.  0.001 .and.
     &         sh(i,j,k) .le.100.001)then
               if(sh(i,j,k).gt.100.)sh(i,j,k)=100.
               sh(i,j,k)=make_ssh(pr(i,j,k)
     .                        ,tp(i,j,k)-273.15
     .                        ,sh(i,j,k)/100.,t_ref)*0.001

               shsum(k)=shsum(k)+sh(i,j,k)
            else
               icm(k)=icm(k)+1
               sh(i,j,k)=rfill
            endif

            if(cmodel(1:nclen).ne.'AVN_FSL_NETCDF'.and.
     +         cmodel(1:nclen).ne.'CWB_20FA_LAMBERT_NF')then
               if(ww(i,j,k).eq.rfill)ww(i,j,k)=0.0
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

c        if(bgmodel.eq.3)then
c           do j=1,ny
c           do i=1,nx
c              mslp(i,j)=mslp(i,j)/100.   !hpa
c           enddo
c           enddo
c        endif

         if(bgmodel.eq.6)then     !no sfc fields for FA model (bgmodel = 3)
 
            print*,'convert rh to Td - sfc: bgmodel: ',bgmodel
            icm_sfc=0
            sumsfctd=0.0
            do j=1,ny
            do i=1,nx

               if(td_sfc(i,j).gt.0.0 .and. td_sfc(i,j).lt.100.001)then
                  prsfc=pr_sfc(i,j)/100.
                  qsfc=make_ssh(prsfc,tp_sfc(i,j)-273.15,td_sfc(i,j)/100.
     &,t_ref)
                  td_sfc(i,j)=make_td(prsfc,tp_sfc(i,j)-273.15,qsfc
     &,t_ref)+273.15
                  sumsfctd=sumsfctd+td_sfc(i,j)
               else
                  td_sfc(i,j)=rfill
c                 td_sfc(i,j)=make_td(pr_sfc(i,j)/100.,tp_sfc(i,j)-273.15
c    &,bogus_sh,t_ref)+273.15
                  icm_sfc=icm_sfc+1
               endif

c           it=tp_sfc(i,j)*100
c           it=min(45000,max(15000,it))
c           xe=esat(it)
c           mrsat=0.00622*xe/(prsfc-xe)         !Assumes that rh units are %
c           td_sfc(i,j)=td_sfc(i,j)*mrsat             !rh --> mr
c           td_sfc(i,j)=td_sfc(i,j)/(1.+td_sfc(i,j))  !mr --> sh

            enddo
            enddo

         endif

      elseif(bgmodel.eq.8)then
c
c *** Convert Td to sh and fill 3D pressure array.
c
         print*,'Convert Td to q - 3D'
         do k=1,nz
          shsum(k)=0.0
          icm(k) = 0
          do j=1,ny
          do i=1,nx
            pr(i,j,k)=prk(k)
            if (sh(i,j,k) .gt. -99999.) then
                sh(i,j,k)=ssh2(prk(k),tp(i,j,k)-273.15,
     & sh(i,j,k)-273.15,t_ref)*0.001

               shsum(k)=shsum(k)+sh(i,j,k)

c              it=sh(i,j,k)*100
c              it=min(45000,max(15000,it))
c              xe=esat(it)
c              sh(i,j,k)=0.622*xe/(pr(i,j,k)-xe)
c              sh(i,j,k)=sh(i,j,k)/(1.+sh(i,j,k))

            else
               sh(i,j,k)=rfill
               icm(k)=icm(k)+1
            endif

            ww(i,j,k)=r_missing_data

          enddo
          enddo
         enddo
c
c check for missing Td NOGAPS data.
c
         icm_sfc=0
         sumsfctd=0.0
         do j=1,ny
         do i=1,nx
            if(td_sfc(i,j).eq.-99999.)then
               icm_sfc=icm_sfc+1
               td_sfc(i,j)=rfill
            else
               sumsfctd=sumsfctd+td_sfc(i,j)
            endif
         enddo
         enddo

      endif
c
c check for and fill missing sh with computed avg for that level
c
      sumtot = 0.0
      icn3d = 0
      do k=1,nz
         if(icm(k).gt.0.and.icm(k).lt.nx*ny)then
            shavg=shsum(k)/(nx*ny-icm(k))
            do j=1,ny
            do i=1,nx
               if(sh(i,j,k).eq.rfill)sh(i,j,k)=shavg
            enddo
            enddo
            sumtot=shsum(k)+sumtot
            icn3d=icm(k)+icn3d
         elseif(icm(k).eq.nx*ny)then
            do j=1,ny
            do i=1,nx
               if(sh(i,j,k).eq.rfill)sh(i,j,k)=r_bogus_sh
            enddo
            enddo
         endif
      enddo

      if(sumtot.gt.0.0.and.icn3d.ne.nx*ny*nz)then
         print*,'Missing background 3D moisture filled with '
     &,'average of good points'
         print*,'#/% 3D points filled: ',icn3d,icn3d/(nx*ny*nz)
      endif

      if(bgmodel.eq.6.or.bgmodel.eq.8)then
         if(icm_sfc.gt.0.and.icm_sfc.lt.nx*ny)then
            shavg=sumsfctd/(nx*ny-icm_sfc)
            do j=1,ny
            do i=1,nx
               if(td_sfc(i,j).eq.rfill)td_sfc(i,j)=shavg
            enddo
            enddo
         endif
      endif
      if(icm_sfc.gt.0)then
         print*,'Missing background 2D moisture filled with '
     &,'average of good points (Td avg = ',shavg,').'
         print*,'#/% 2D points filled: ',icm_sfc,icm_sfc/(nx*ny)
      endif

      istatus=0
      return
c

990   print *,'Error finding readdgprep file.'

      IF (IOSTATUS .NE. 0)THEN
         PRINT *,'ERROR READING ',FILENAME(1:l),' IO status is',
     &  IOSTATUS
      END IF

      inquire(file=filename,exist=lext,opened=lopen,number=n)
      if(.not.lext)then
         print*,'File does not exist: ',filename(1:l)
      endif
      if(lopen)then
         print*,'File is already open: ',filename(1:l)
      endif
      if(n .lt. 0)then
         print*,'No unit number associated with file '
     &                ,filename(1:l)
      endif

      return
995   print*,'Error reading model index file.',filename(1:l)
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
c qc for model rh=0.0. 
c
      call qcmodel_sh(nx,ny,1,sh_sfc)

      call qcmodel_sh(nx,ny,nz,sh)

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
