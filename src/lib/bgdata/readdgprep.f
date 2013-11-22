      subroutine read_dgprep(bgmodel,cmodel,path,fname,af,nx,ny,nz
     .                      ,pr,ht,tp,sh,uw,vw,ww
     .                      ,ht_sfc,pr_sfc,sh_sfc,tp_sfc,t_at_sfc
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
     .         ,icn3d,ipk(nz)
c
      real     shsum(nz),sumtot,sumsfctd
      real     shavg
      integer  it
      integer  lun
      integer  iostat,iostatus
      integer  nclen,nflen,lenc
      integer  version

      logical  lopen,lext,ldo_tcbogus

      real   ht(nx,ny,nz)      !height (m)
     .      ,tp(nx,ny,nz)      !temperature (K)
     .      ,sh(nx,ny,nz)      !specific humidity (kg/kg) 
     .      ,uw(nx,ny,nz)      !u-wind (m/s)
     .      ,vw(nx,ny,nz)      !v-wind (m/s)
     .      ,ww(nx,ny,nz)      !w-wind (pa/s)  !currently not available for NOGAPS
     .      ,pr(nx,ny,nz)      !pressures (mb)
     .      ,pw(nx,ny,nz)      !precip h2o for /public AVN
     .      ,prk(nz)

      real   ht_sfc(nx,ny)
     .      ,pr_sfc(nx,ny)     !in pascals
     .      ,sh_sfc(nx,ny)     !specific humidity (kg/kg) 
     .      ,tp_sfc(nx,ny)     !Deg K, for /public AVN this is T @ 2m agl.
     .      ,rh_sfc(nx,ny)     !Percent
     .      ,uw_sfc(nx,ny)
     .      ,vw_sfc(nx,ny)
     .      ,mslp(nx,ny)
     .      ,accs(nx,ny)       !accum snow for /public AVN
     .      ,t_at_sfc(nx,ny)   !only used for /public AVN  and CWB NF15... holds T @ sfc
        
      real   p_levels(nz,nvarsmax)

      double precision isoLevel(nz),reftime,valtime

c     real   mrsat
c     real   esat,xe
c     real   rp_init

      real   prsfc,prsfc_mb
      real   qsfc
      real   make_td
      real   make_ssh
      real   r_bogus_sh
      real   ssh2
      real   t_ref
      real   pcnt
      real   r_missing_data
      real   rfill
      real   rmx2d,rmn2d

      real   lat1,lat2,lon0,dskm,glat(nx,ny),glon(nx,ny),sw(2),ne(2)       

      integer imx,jmx,imn,jmn

c
      character*256 path
      character*132 cmodel
      character*200 fname
      character*10  cfname10,a9_to_yr_a10_time
      character*4   af
      character*2   gproj
      character*2   cpref
 
      character*16  cFA_filename
      character*3   c3ext,  c3_FA_ext
      character*8   cwb_type
      character*132 origin,model,nav,grid
      character*255 filename,fname_index
c
c     common /estab/esat(15000:45000)
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
      istatus = 0
      icm_sfc = 0
c
c This is set = .true. only for Taiwan runs ATM.
c 
      ldo_tcbogus=.false.
c
      call get_r_missing_data(r_missing_data,istatus)

!     Initialize surface fields
      ht_sfc = r_missing_data
      pr_sfc = r_missing_data
      tp_sfc = r_missing_data
      sh_sfc = r_missing_data
      rh_sfc = r_missing_data
      uw_sfc = r_missing_data
      vw_sfc = r_missing_data

      call s_len(cmodel,nclen)
      call s_len(fname,nflen)

      lun=10
c
      if(bgmodel.eq.6.or.bgmodel.eq.8.or.bgmodel.eq.12)then

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
     +,ht_sfc,pr_sfc,sh_sfc,tp_sfc,uw_sfc,vw_sfc,mslp
     +,istatus)

           elseif(bgmodel.eq.8)then

              call read_nogaps(lun,nx,ny,nz
     + ,nvarsmax,nvars,nlevs,ivarcoord,ivarid
     + ,ht,tp,sh,uw,vw,ht_sfc,pr_sfc,sh_sfc,tp_sfc
     + ,uw_sfc,vw_sfc,mslp,istatus)


           endif

c        else
c
c eta ingest currently disabled. J. Smart (9-2-98)
c
c           call read_eta(lun,nx,ny,nz,tp,uw,vw,ht,sh
c    +,ht_sfc,pr_sfc,sh_sfc,tp_sfc,uw_sfc,vw_sfc,mslp
c    +,istatus)

         elseif(cmodel(1:nclen).eq.'AVN_FSL_NETCDF')then

               call read_avn_netcdf(filename, nz, 1, nx, ny,
     +     version, ACCS, ht , ht_sfc, pr_sfc, mslp, pw, sh, tp,
     +     tp_sfc, t_at_sfc, uw, vw, ww, sh_sfc, isoLevel,
     +     reftime, valtime, grid, model, nav, origin, istatus)


               nzsh=nz-5

               call qcmodel_sh(nx,ny,1,sh_sfc)  !sh_sfc actually = RH for AVN.

               call qcmodel_sh(nx,ny,nzsh,sh)     !sh actually = RH for AVN.

         elseif(cmodel(1:nclen).eq.'FMI_NETCDF_LL')then

              call read_FMI_netcdf(filename, nz, 1, nx, ny,
     +     ipk,sh_sfc,ht,ht_sfc,mslp,pr_sfc,ww,sh,tp,tp_sfc,t_at_sfc,
     +     uw, uw_sfc, vw, vw_sfc, model, origin, reftime, valtime,
     +     istatus)
              if(istatus.ne.0)then
               print*,'Error!!! Model read failure!!!'
               print*,'Filename ',TRIM(filename)
               return
              endif
              prk=float(ipk)
              nzsh=nz
         endif

      elseif(bgmodel.eq.3)then

         ldo_tcbogus=.true.
c
c note: library function fname13_to_FA_filename could be used
c       to covert the FA filename but currently is not.  J.Smart

         search_mdltype: do l=nclen,1,-1
            if(cmodel(l:l).eq.'_')then
               exit search_mdltype
            endif
         enddo search_mdltype
         cwb_type=cmodel(l+1:nclen)

         call downcase(cwb_type,cwb_type)
         c3ext=c3_FA_ext(af) 
         cfname10=a9_to_yr_a10_time(fname,istatus)
         if(cwb_type.eq.'gfs')then
            cpref='gb'
         elseif(cwb_type.eq.'tfs')then
            cpref='sb'
         else
            cpref='nf'
         endif
cFA_filename="nf"//cfname10(1:8)//fname(6:7)//'.'//c3ext
         cFA_filename=cpref//cfname10(1:8)//fname(6:7)//'.'//c3ext
         call s_len(path,l)
         filename=path(1:l)//'/'//cFA_filename
         call s_len(filename,l)

         if(cwb_type .eq. 're')then
            print*,'Opening FA file: ',filename(1:l)
            open(lun,file=filename(1:l),status='old'
     +,IOSTAT=IOSTATUS,err=990)

            call read_fa(lun,filename                   ! I
     .               ,nx,ny,nz                          ! I
     .               ,r_missing_data                    ! I
     .               ,prk                               ! O
     .               ,ht,tp,sh,uw,vw                    ! O
     .               ,mslp                              ! O
     .               ,istatus)                          ! O

            close(lun)
            call qcmodel_sh(nx,ny,nz,sh)
         elseif (cwb_type.eq.'nf')then
            print *,'read nf model with read_fa_nf' 
            open(lun,file=filename(1:l),status='old'
     +,IOSTAT=IOSTATUS,err=990)

            call read_fa_nf(lun,filename,              ! I
     .                nx, ny, nz,                      ! I
     .                r_missing_data,                  ! I
     .                prk,                             ! O
     .                ht, tp, sh, uw, vw, ww,          ! O
     .                ht_sfc, tp_sfc, rh_sfc,          ! O
     .                uw_sfc, vw_sfc, mslp,            ! O
     .                istatus )                        ! O
            close(lun)
            call qcmodel_sh(nx,ny,nz,sh)
c           do k=1,nz
c              call get_mxmn_2d(nx,ny,ww,rmx2d,rmn2d
c    &,imx,jmx,imn,jmn)
c              print*,'k mx/mn ww ',k,rmx2d,rmn2d
c           enddo

            ww=ww/36.
            pr_sfc = pr_sfc * 100. ! convert to Pascals

         elseif (cwb_type.eq.'nf15'.or.cwb_type.eq.'gfs'
     &.or.cwb_type.eq.'nf45'.or.cwb_type.eq.'tfs')then
            print *,'read nf model with read_nf15km' 

            t_at_sfc = r_missing_data

            call read_nf15km(nx,ny,nz,filename,
     &                       ht,tp,sh,uw,vw,ww,        !Note: sh contains 3D rh
     &                       pr_sfc,tp_sfc,rh_sfc,    
     &                       uw_sfc,vw_sfc,mslp,
     &                       prk,t_at_sfc,
     &                       istatus)

            ww=ww/36.              ! convert to pa/s
            pr_sfc = pr_sfc * 100. ! convert to Pascals

         else
            print*,'cwb model type currently unknown ',cmodel(1:nclen)
            print*,'return without data'
            return
         endif

         nzsh=nz

      endif

      if(istatus .eq. 1)then
         print*,'Error reading dgprep data file: ',cmodel(1:nclen),
     +' ',filename(1:l)
         return
      endif

      if(ldo_tcbogus)then

!        Add call to get_bkgd_mdl_info for lat1,lat2,lon0,sw,ne,dskm

!        Add call to init_hinterp to get glat,glon

!        Right now the bogusing assumes a Lambert Conformal Projection

         call tcbogus(nx,ny,nz,ht,tp,sh,uw,vw,ht_sfc,
     +             tp_sfc,sh_sfc,uw_sfc,vw_sfc,mslp,
     +             lat1,lat2,lon0,sw,ne,dskm,glat,glon,
     +             prk,filename,bgmodel,cwb_type)

      endif

c ---------------------------------------------
c qc 
c --
      do k=1,nz
         if(k.le.nzsh)then
            ksh=k
         else
            ksh=nzsh
         endif
         do j=1,ny
         do i=1,nx
            if((abs(ht(i,j,k))  .gt. 110000.) .or.
     +             (ht(i,j,k)   .lt.-3000.)   .or.
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
c *** Note: sh and sh_sfc arrays contain rh from AVN (bgmodel=6)
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

!              Check if sh_sfc contains valid RH value
               if(sh_sfc(i,j).gt.0.0 .and. sh_sfc(i,j).lt.100.001)then
                  prsfc=pr_sfc(i,j)/100.
                  qsfc=make_ssh(prsfc,tp_sfc(i,j)-273.15,sh_sfc(i,j)/100.
     &,t_ref)
                  sh_sfc(i,j)=make_td(prsfc,tp_sfc(i,j)-273.15,qsfc
     &,t_ref)+273.15
                  sumsfctd=sumsfctd+sh_sfc(i,j)
               else 
                  sh_sfc(i,j)=rfill
c                 sh_sfc(i,j)=make_td(pr_sfc(i,j)/100.,tp_sfc(i,j)-273.15
c    &,bogus_sh,t_ref)+273.15
                  icm_sfc=icm_sfc+1
               endif

c           it=tp_sfc(i,j)*100
c           it=min(45000,max(15000,it))
c           xe=esat(it)
c           mrsat=0.00622*xe/(prsfc-xe)         !Assumes that rh units are %
c           sh_sfc(i,j)=sh_sfc(i,j)*mrsat             !rh --> mr
c           sh_sfc(i,j)=sh_sfc(i,j)/(1.+sh_sfc(i,j))  !mr --> sh

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
            if(sh_sfc(i,j).eq.-99999.)then
               icm_sfc=icm_sfc+1
               sh_sfc(i,j)=rfill
            else
               sumsfctd=sumsfctd+sh_sfc(i,j)
            endif
         enddo
         enddo

      elseif(bgmodel.eq.12.and.
     &cmodel(1:nclen).eq.'FMI_NETCDF_LL')then
        do k=1,nz
          do j=1,ny
           do i=1,nx
            pr(i,j,k)=prk(k)
           enddo
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
               if(sh_sfc(i,j).eq.rfill)sh_sfc(i,j)=shavg
            enddo
            enddo
         endif
      endif
      if(icm_sfc.gt.0)then
         print*,'Missing background 2D moisture filled with '
     &,'average of good points (Td avg = ',shavg,').'
         print*,'#/% 2D points filled: ',icm_sfc,icm_sfc/(nx*ny)
      endif

!     Convert rh_sfc to sh_sfc if sh_sfc is unavailable
      if(minval(sh_sfc) .eq. r_missing_data .OR. 
     1   maxval(sh_sfc) .eq. r_missing_data      )then
          write(6,*)' Calculate sh_sfc using Sfc T, P, RH'
          do i = 1,nx
          do j = 1,ny
              if(rh_sfc(i,j) .ne. r_missing_data .AND.
     1           pr_sfc(i,j) .ne. r_missing_data .AND.
     1           tp_sfc(i,j) .ne. r_missing_data .AND.
     1           rh_sfc(i,j) .ge. 0.0            .AND.
     1           rh_sfc(i,j) .le. 100.0                )then
                  prsfc_mb=pr_sfc(i,j)/100.
                  sh_sfc(i,j)=(make_ssh(prsfc_mb,tp_sfc(i,j)-273.15
     1                        ,rh_sfc(i,j)/100.,t_ref)) * .001
!                 td_sfc(i,j)=make_td(prsfc_mb,tp_sfc(i,j)-273.15,qsfc
!    1                               ,t_ref)+273.15
              endif
          enddo ! j
          enddo ! i
      endif

      write(6,*)' Returning from read_dgprep'

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
c
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

      real   ht(nx,ny,nz)      !height (m)
     .      ,tp(nx,ny,nz)      !temperature (K)
     .      ,sh(nx,ny,nz)      !specific humidity (kg/kg)
     .      ,uw(nx,ny,nz)      !u-wind (m/s)
     .      ,vw(nx,ny,nz)      !v-wind (m/s)

      real   ht_sfc(nx,ny)
     .      ,pr_sfc(nx,ny)
     .      ,sh_sfc(nx,ny)
     .      ,tp_sfc(nx,ny)
     .      ,uw_sfc(nx,ny)
     .      ,vw_sfc(nx,ny)
     .      ,mslp(nx,ny)

      real   dummy(nx,ny,nz)

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
c
      subroutine read_eta(lun,nx,ny,nz,tp,uw,vw,ht,sh
     +,ht_sfc,pr_sfc,sh_sfc,tp_sfc,uw_sfc,vw_sfc,mslp
     +,istatus)

      implicit none

      integer   nx,ny,nz
     .         ,i,j,k,istatus
c
      integer  lun
      real   ht(nx,ny,nz)      !height (m)
     .      ,tp(nx,ny,nz)      !temperature (K)
     .      ,sh(nx,ny,nz)      !specific humidity (kg/kg)
     .      ,uw(nx,ny,nz)      !u-wind (m/s)
     .      ,vw(nx,ny,nz)      !v-wind (m/s)

      real   ht_sfc(nx,ny)
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
C ---------------------------------------------------------------
C
      subroutine qcmodel_sh(nx,ny,nz,sh)

      implicit none

      integer i,j,k,nx,ny,nz

      real    sh(nx,ny,nz)

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

! ---------------------------------------------------------

      subroutine read_nf15km(mx,my,nz,full_name
     &,ht_ou,tp_ou,rh_ou,uw_ou,vw_ou,ww_ou
     &,pss_ou,tps_ou,rhs_ou,uws_ou,vws_ou,mslp_ou
     &,prk,tmp_ou,istatus)

c
c  GFS and NEW NFS_15km model information 
c  
c  GFS had been interpolatted to NEW NFS_15 same domain and projection
c
c               mx=181,my=193,nz=11 
c  SW corner (9.28194 N, 109.7727E), NE corner (35.26665N,137.7342E), 
c  map projection same as now NF model . 
c    nz : 1-> 1000mb,2-> 925mb ,3-> 850mb,4->700mb, 5 -> 500mb, 6-> 400mb
c    nz : 7-> 300mb,8-> 250mb ,9-> 200mb,10->150mb,11 -> 100mb
c 
c
c               3D Output fields (pressure grid)
      REAL      :: ht_ou(mx,my,nz)     ! nfs height (m)
      REAL      :: tp_ou(mx,my,nz)     ! nfs temperature (K)
      REAL      :: rh_ou(mx,my,nz)     ! nfs relative humidity (%)
      REAL      :: uw_ou(mx,my,nz)     ! nfs u-wind (m/s)
      REAL      :: vw_ou(mx,my,nz)     ! nfs v-wind (m/s)
      REAL      :: ww_ou(mx,my,nz)     ! nfs w-wind (hPa/hr)
c!                    2D Output fields (Sfc field)
      REAL      :: pss_ou(mx,my)       ! nfs surface pressure (hPa)
      REAL      :: tps_ou(mx,my)       ! nfs temperature (K)
      REAL      :: rhs_ou(mx,my)       ! nfs relative humidity (%)
      REAL      :: uws_ou(mx,my)       ! nfs u-wind (m/s)
      REAL      :: vws_ou(mx,my)       ! nfs v-wind (m/s)
      REAL      :: mslp_ou(mx,my)      ! nfs mean sea level presures (pa)
      REAL      :: tmp_smt(mx,my)
c                     pressures of the levels
      REAL      :: prk(nz) 	       ! pressure of each level (pa)
c
      CHARACTER(LEN=256)  :: header_ht
      CHARACTER(LEN=256)  :: header_tp
      CHARACTER(LEN=256)  :: header_rh
      CHARACTER(LEN=256)  :: header_uw
      CHARACTER(LEN=256)  :: header_vw
      CHARACTER(LEN=256)  :: header_ww
      CHARACTER(LEN=256)  :: header_pss
      CHARACTER(LEN=256)  :: header_tps
      CHARACTER(LEN=256)  :: header_rhs
      CHARACTER(LEN=256)  :: header_uws
      CHARACTER(LEN=256)  :: header_vws
      CHARACTER(LEN=256)  :: header_mslp

      CHARACTER*(*)       :: full_name

      integer istatus,l
c
     
c
      header_ht  ='nfs height (m)'
      header_tp  ='nfs temperature (K)'
      header_rh  ='nfs relative humidity (%)'
      header_uw  ='nfs u-wind (m/s)'
      header_vw  ='nfs v-wind (m/s)'
      header_ww  ='nfs w-wind (hPa/hr)'
      header_pss ='nfs surface pressure (hPa)'
      header_tps ='nfs temperature (K)'
      header_rhs ='nfs relative humidity (%)'
      header_uws ='nfs u-wind (m/s)'
      header_vws ='nfs v-wind (m/s)'
      header_mslp='nfs mean sea level presures (hpa)'  !modified to "hpa" 9-22-04 

      prk(1) = 100000.
      prk(2) = 92500.
      prk(3) = 85000.
      prk(4) = 70000.
      prk(5) = 50000.
      prk(6) = 40000.
      prk(7) = 30000.
      prk(8) = 25000.
      prk(9) = 20000.
      prk(10)= 15000.
      prk(11)= 10000.

c   header and  data
c
      call s_len(full_name,l)
      istatus = 1
      open(18,file=full_name,status='old',form='unformatted')
      read(18,err=10)header_ht
      read(18,err=10)ht_ou
      read(18,err=10)header_tp
      read(18,err=10)tp_ou
      read(18,err=10)header_rh
      read(18,err=10)rh_ou
      read(18,err=10)header_uw
      read(18,err=10)uw_ou
      read(18,err=10)header_vw
      read(18,err=10)vw_ou
      read(18,err=10)header_ww
      read(18,err=10)ww_ou
      read(18,err=10)header_pss
      read(18,err=10)pss_ou
      read(18,err=10)header_tps
      read(18,err=10)tps_ou
      read(18,err=10)header_rhs
      read(18,err=10)rhs_ou
      read(18,err=10)header_uws
      read(18,err=10)uws_ou
      read(18,err=10)header_vws
      read(18,err=10)vws_ou
      read(18,err=10)header_mslp
      read(18,err=10)mslp_ou
      close(18)

      istatus = 0
      return

10    print*,'Error reading file ',full_name(1:l)
      return
      end


