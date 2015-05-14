
      subroutine terrain_scale(nx,ny,topo,terscl)

      implicit none
      integer  nx,ny

      integer  imx,jmx,imn,jmn
      integer  i,j

      real     topo(nx,ny)
      real     rmx2d,rmn2d
      real     sumt
      real     terscl

      call get_mxmn_2d(nx,ny,topo,rmx2d,rmn2d
     &,imx,jmx,imn,jmn)

      terscl=rmx2d-rmn2d

      sumt=0.
      do j=1,ny
       do i=1,nx
        sumt=sumt+( (topo(i,j)-rmn2d)**2
     &+(topo(i,j)-rmx2d)**2)
       enddo
      enddo
      terscl=sqrt(sumt/(2.*nx*ny))

      return
      end
c ----------------------------------------------------------
c more sophisticated routine for variable tau and k.
c
      subroutine terrain_scale_vartau(nx,ny,nz,topo,hgt3d
     &,ks,kf,terscl)

      implicit none
      integer  nx,ny,nz

c     integer  imx,jmx,imn,jmn
      integer  i,j,k,ii,jj
      integer  ks(nx,ny)
      integer  kf(nx,ny)
      integer  istatus

      logical  lfhf,lfhs

      real     topo(nx,ny)
      real     hgt3d(nx,ny,nz)
      real     rmxt,rmnt
      real     r_missing_data
c     real     sumt
c     real     dif1,dif2
c     real     toposcl

      real     terscl(nx,ny)

c     call get_mxmn_2d(nx,ny,topo,rmx2d,rmn2d
c    &,imx,jmx,imn,jmn)

c     terscl=rmx2d-rmn2d

c     sumt=0.

      call get_r_missing_data(r_missing_data,istatus)

      do j=2,ny-1
       do i=2,nx-1

c       dif1=topo(i,j)-rmn2d
c       dif2=topo(i,j)-rmx2d
c       sumt=sumt+(dif1*dif1)+(dif2*dif2) 

        rmxt=-1*r_missing_data
        rmnt=r_missing_data
        do jj=j-1,j+1
        do ii=i-1,i+1
         rmxt=max(rmxt,topo(ii,jj))
         rmnt=min(rmnt,topo(ii,jj))
        enddo
        enddo
        terscl(i,j)=rmxt-rmnt
        if(terscl(i,j).lt.1.0)terscl(i,j)=1.0
        
        lfhf=.false.
        lfhs=.false.
        do k=1,nz
         if((hgt3d(i,j,k).gt.rmxt).and.(.not.lfhf))then
          kf(i,j)=k+1
          lfhf=.true.
         endif
         if((hgt3d(i,j,k).gt.rmnt).and.(.not.lfhs))then
          ks(i,j)=k
          if(k.eq.1)ks(i,j)=k
          lfhs=.true.
         endif
        enddo

       enddo
      enddo

c N/S boundaries
      do i=2,nx-1
       rmxt=-1*r_missing_data
       rmnt=r_missing_data
       do j=1,2
        rmxt=max(rmxt,topo(i-1,j))
        rmnt=min(rmnt,topo(i-1,j))
        rmxt=max(rmxt,topo(i,j))
        rmnt=min(rmnt,topo(i,j))
        rmxt=max(rmxt,topo(i+1,j))
        rmnt=min(rmnt,topo(i+1,j))
       enddo
       terscl(i,1)=rmxt-rmnt
       if(terscl(i,1).lt.1.0)terscl(i,1)=1.0
c
c k values
c
       do j=1,2
        lfhf=.false.
        lfhs=.false.
        do k=1,nz
         if(((hgt3d(i-1,j,k).gt.rmxt).or.
     &       (hgt3d(i,j,k)  .gt.rmxt).or.
     &       (hgt3d(i+1,j,k).gt.rmxt)).and.(.not.lfhf))then
             kf(i,j)=k+1
             lfhf=.true.
         endif
         if(((hgt3d(i-1,1,k).gt.rmnt).or.
     &       (hgt3d(i,1,k)  .gt.rmnt).or.
     &       (hgt3d(i+1,1,k).gt.rmnt)).and.(.not.lfhs))then
             ks(i,1)=k
             if(k.eq.1)ks(i,1)=k
             lfhs=.true.
         endif
        enddo
       enddo

       rmxt=-1*r_missing_data
       rmnt=r_missing_data
       do j=ny-1,ny
        rmxt=max(rmxt,topo(i-1,j))
        rmnt=min(rmnt,topo(i-1,j))
        rmxt=max(rmxt,topo(i,j))
        rmnt=min(rmnt,topo(i,j))
        rmxt=max(rmxt,topo(i+1,j))
        rmnt=min(rmnt,topo(i+1,j))
       enddo
       terscl(i,ny)=rmxt-rmnt
       if(terscl(i,ny).lt.1.0)terscl(i,ny)=1.0 
c
c k values
c
       do j=ny-1,ny
        lfhf=.false.
        lfhs=.false.
        do k=1,nz
         if(((hgt3d(i-1,j,k).gt.rmxt).or.
     &       (hgt3d(i,j,k)  .gt.rmxt).or.
     &       (hgt3d(i+1,j,k).gt.rmxt)).and.(.not.lfhf))then
             kf(i,ny)=k+1
             lfhf=.true.
         endif
         if(((hgt3d(i-1,j,k).gt.rmnt).or.
     &       (hgt3d(i,j,k)  .gt.rmnt).or.
     &       (hgt3d(i+1,j,k).gt.rmnt)).and.(.not.lfhs))then
             ks(i,ny)=k
             if(k.eq.1)ks(i,ny)=k
             lfhs=.true.
         endif
        enddo
       enddo

      enddo
c
c E/W boundaries
      do j=2,ny-1
       rmxt=-1*r_missing_data
       rmnt=r_missing_data
       do i=1,2
        rmxt=max(rmxt,topo(i,j-1))
        rmnt=min(rmnt,topo(i,j-1))
        rmxt=max(rmxt,topo(i,j))
        rmnt=min(rmnt,topo(i,j))
        rmxt=max(rmxt,topo(i,j+1))
        rmnt=min(rmnt,topo(i,j+1))
       enddo
       terscl(1,j)=rmxt-rmnt
       if(terscl(1,j).lt.1.0)terscl(1,j)=1.0
c
c k values
c
       do i=1,2
        lfhf=.false.
        lfhs=.false.
        do k=1,nz
         if(((hgt3d(i,j-1,k).gt.rmxt).or.
     &       (hgt3d(i,j,k)  .gt.rmxt).or.
     &       (hgt3d(i,j+1,k).gt.rmxt)).and.(.not.lfhf))then
             kf(1,j)=k+1
             lfhf=.true.
         endif
         if(((hgt3d(i,j-1,k).gt.rmnt).or.
     &       (hgt3d(i,j,k)  .gt.rmnt).or.
     &       (hgt3d(i,j+1,k).gt.rmnt)).and.(.not.lfhs))then
             ks(1,j)=k
             if(k.eq.1)ks(1,j)=k
             lfhs=.true.
         endif
        enddo
       enddo

       rmxt=-1*r_missing_data
       rmnt=r_missing_data
       do i=nx-1,nx
        rmxt=max(rmxt,topo(i,j-1))
        rmnt=min(rmnt,topo(i,j-1))
        rmxt=max(rmxt,topo(i,j))
        rmnt=min(rmnt,topo(i,j))
        rmxt=max(rmxt,topo(i,j+1))
        rmnt=min(rmnt,topo(i,j+1))
       enddo
       terscl(nx,j)=rmxt-rmnt
       if(terscl(nx,j).lt.1.0)terscl(nx,j)=1.0
c
c k values
c
       do i=nx-1,nx
        lfhf=.false.
        lfhs=.false.
        do k=1,nz
         if(((hgt3d(i,j-1,k).gt.rmxt).or.
     &       (hgt3d(i,j,k)  .gt.rmxt).or.
     &       (hgt3d(i,j+1,k).gt.rmxt)).and.(.not.lfhf))then
             kf(nx,j)=k+1
             lfhf=.true.
         endif
         if(((hgt3d(i,j-1,k).gt.rmnt).or.
     &       (hgt3d(i,j,k)  .gt.rmnt).or.
     &       (hgt3d(i,j+1,k).gt.rmnt)).and.(.not.lfhs))then
             ks(nx,j)=k
             if(k.eq.1)ks(nx,j)=k
             lfhs=.true.
         endif
        enddo
       enddo

      enddo

c corners
c SW/NW
      rmxt=-1*r_missing_data
      rmnt=r_missing_data
      do i=1,2
       do j=1,2
        rmxt=max(rmxt,topo(i,j))
        rmnt=min(rmnt,topo(i,j))
       enddo
      enddo
      terscl(1,1)=rmxt-rmnt
      if(terscl(1,1).lt.1.0)terscl(1,1)=1.0
c
c k values
c
      do i=1,2
        lfhf=.false.
        lfhs=.false.
        do k=1,nz
         if(((hgt3d(i,1,k).gt.rmxt).or.
     &       (hgt3d(i,2,k).gt.rmxt)).and.(.not.lfhf))then
             kf(1,1)=k+1
             lfhf=.true.
         endif
         if(((hgt3d(i,1,k).gt.rmnt).or.
     &       (hgt3d(i,2,k).gt.rmnt)).and.(.not.lfhs))then
             ks(1,1)=k
             if(k.eq.1)ks(1,1)=k
             lfhs=.true.
         endif
        enddo
      enddo

      rmxt=-1*r_missing_data
      rmnt=r_missing_data
      do i=1,2
       do j=ny-1,ny
        rmxt=max(rmxt,topo(i,j))
        rmnt=min(rmnt,topo(i,j))
       enddo
      enddo
      terscl(1,ny)=rmxt-rmnt
      if(terscl(1,ny).lt.1.0)terscl(1,ny)=1.0
c
c k values
c
      do i=1,2
        lfhf=.false.
        lfhs=.false.
        do k=1,nz
         if(((hgt3d(i,ny,k)  .gt.rmxt).or.
     &       (hgt3d(i,ny-1,k).gt.rmxt)).and.(.not.lfhf))then
             kf(1,ny)=k+1
             lfhf=.true.
         endif
         if(((hgt3d(i,ny,k)  .gt.rmnt).or.
     &       (hgt3d(i,ny-1,k).gt.rmnt)).and.(.not.lfhs))then
             ks(1,ny)=k
             if(k.eq.1)ks(1,ny)=k
             lfhs=.true.
         endif
        enddo
      enddo

c SE/NE
      rmxt=-1*r_missing_data
      rmnt=r_missing_data
      do i=nx-1,nx
       do j=1,2
        rmxt=max(rmxt,topo(i,j))
        rmnt=min(rmnt,topo(i,j))
       enddo
      enddo
      terscl(nx,1)=rmxt-rmnt
      if(terscl(nx,1).lt.1.0)terscl(nx,1)=1.0
c
c k values
c
      do i=nx-1,nx
        lfhf=.false.
        lfhs=.false.
        do k=1,nz
         if(((hgt3d(i,1,k).gt.rmxt).or.
     &       (hgt3d(i,2,k).gt.rmxt)).and.(.not.lfhf))then
             kf(nx,1)=k+1
             lfhf=.true.
         endif
         if(((hgt3d(i,1,k).gt.rmnt).or.
     &       (hgt3d(i,2,k).gt.rmnt)).and.(.not.lfhs))then
             ks(nx,1)=k
             if(k.eq.1)ks(nx,1)=k
             lfhs=.true.
         endif
        enddo
      enddo

      rmxt=-1*r_missing_data
      rmnt=r_missing_data
      do i=nx-1,nx
       do j=ny-1,ny
        rmxt=max(rmxt,topo(i,j))
        rmnt=min(rmnt,topo(i,j))
       enddo
      enddo
      terscl(nx,ny)=rmxt-rmnt
      if(terscl(nx,ny).lt.1.0)terscl(nx,ny)=1.0
c
c k values
c
      do i=nx-1,nx
        lfhf=.false.
        lfhs=.false.
        do k=1,nz
         if(((hgt3d(i,ny,k).gt.rmxt).or.
     &       (hgt3d(i,ny-1,k).gt.rmxt)).and.(.not.lfhf))then
             kf(nx,ny)=k+1
             lfhf=.true.
         endif
         if(((hgt3d(i,ny,k).gt.rmnt).or.
     &       (hgt3d(i,ny-1,k).gt.rmnt)).and.(.not.lfhs))then
             ks(nx,ny)=k
             if(k.eq.1)ks(nx,ny)=k
             lfhs=.true.
         endif
        enddo
      enddo

c     terscl=sqrt(sumt/(nx*ny))

      return
      end

      subroutine advance_grids(i4time_sys_sonde,i4time_sys_drop
     .,nx,ny,nz,ub_a,vb_a,tb_a,phib_a,shb_a,omb_a,psb_a
     .,lapsphi,lapstemp,lapsu,lapsv,lapssh,lapsomo,ps,istatus)

      implicit none

      include 'bgdata.inc'

      integer nx,ny,nz

      real    phib_a(nx,ny,nz)
     .       ,tb_a(nx,ny,nz)
     .       ,ub_a(nx,ny,nz)
     .       ,vb_a(nx,ny,nz)
     .       ,shb_a(nx,ny,nz)
     .       ,omb_a(nx,ny,nz)
     .       ,psb_a(nx,ny)

      real    lapsu(nx,ny,nz)
     .       ,lapsv(nx,ny,nz)
     .       ,lapssh(nx,ny,nz)
     .       ,lapstemp(nx,ny,nz)
     .       ,lapsphi(nx,ny,nz)
     .       ,lapsomo(nx,ny,nz)
     .       ,ps(nx,ny)

      real, allocatable, dimension(:,:,:) ::
     .      phib, tb, ub, vb, shb, omb
      real, allocatable :: psb(:,:)
 

      integer   bgmodels(maxbgmodels)
      integer   itime_inc,ntmin,ntmax

      character(len=256), allocatable, dimension(:) ::
     .        names,reject_names,bg_names,bgpaths
      character(len=132), allocatable, dimension(:) :: cmodels

      integer   max_files
      parameter (max_files=200)
      integer accepted_files
      integer bg_files
      integer forecast_length
      integer rejected_cnt
      integer laps_cycle_time
      logical use_analysis
      logical use_forecast
      logical luse_sfc_bkgd
      logical smooth_fields
      logical lgb_only

      integer i4time_sys_sonde
     .       ,i4time_sys_drop
     .       ,lenm
     .       ,istatus
     .       ,lga_status

c-------------------------------------------------------
c                 \
c              ^=====>  (<- airplane)
c                 /
c
c                dropsonde          payload drop
c                  (ta)                (td)
c                    ^
c                    !                   !
c            bkgd    !   bkgd    bkgd    !   bkgd
c              1     !     2       3     !     4
c        - - - - - - - - - - - - - - - - - - - - - 
c              |     !     |       |frac1!frac2|
c                    !                   !
c                    ^                   ^
c                 analysis            airdrop
c                   (A)                 (D)
c                 input:lga,anal
c 
c ----------------------- time -> ---------------------
c    background at D (this happens in lga_driver):
c       bkgd_lga-D = (frac2*bkgd4 + frac1*bkgd3)
c
c        where bkgd is defined in background.nl (lga namelist)
c
c    analysis at D:
c       analysis-D = analysis-A + (bkgd_lga-D - input_lga-A)
c
c    return bkgd_lga-D and analysis-D for qbalpe.
c
      print*,'  Subroutine advance_grids'

      call get_laps_cycle_time(laps_cycle_time,istatus)
      if(istatus.ne.1)then
         print*,'error returned: no cycle time. return'
         return
      endif

      allocate (cmodels(maxbgmodels),bgpaths(maxbgmodels))

      call get_background_info(bgpaths,bgmodels,forecast_length
     +,use_analysis,use_forecast,cmodels,itime_inc,smooth_fields
     +,luse_sfc_bkgd
     +,ntmin,ntmax
     +,lgb_only)

      call s_len(cmodels(1),lenm)

      allocate (names(max_files),reject_names(max_files)
     .         ,bg_names(max_files))

      print*,'Calling get_acceptable_files'

      call get_acceptable_files(i4time_sys_drop,bgpaths(1),bgmodels(1)
     +        ,names,max_files
     +        ,use_analysis,use_forecast,bg_files,accepted_files
     +        ,forecast_length
     +        ,cmodels(1),nx,ny,nz,reject_names,rejected_cnt)

      print*,'Returned from get_acceptable_files'

      print*,'Calling lga_driver!'
      call lga_driver(nx,ny,nz,luse_sfc_bkgd,laps_cycle_time
     .         ,bgmodels(1),bgpaths(1),cmodels(1),rejected_cnt
     .         ,reject_names,names,max_files,accepted_files
     .         ,i4time_sys_drop,smooth_fields,lgb_only,lga_status)

      if(lga_status.ne.1)then
         print*,'error returned from lga_driver'
         print*,'terminating in advance_grids'
         istatus = 0
         return
      endif

      deallocate (names,reject_names,bg_names,cmodels,bgpaths)
c
c now read the lga/lgb just made
c
      allocate (ub(nx,ny,nz)
     .         ,vb(nx,ny,nz)
     .         ,tb(nx,ny,nz)
     .         ,phib(nx,ny,nz)
     .         ,shb(nx,ny,nz)
     .         ,omb(nx,ny,nz)
     .         ,psb(nx,ny))

      call get_modelfg_3d(i4time_sys_drop,'U3 ',nx,ny,nz,ub,istatus)
      call get_modelfg_3d(i4time_sys_drop,'V3 ',nx,ny,nz,vb,istatus)
      call get_modelfg_3d(i4time_sys_drop,'T3 ',nx,ny,nz,tb,istatus)
      call get_modelfg_3d(i4time_sys_drop,'HT ',nx,ny,nz,phib,istatus)
      call get_modelfg_3d(i4time_sys_drop,'SH ',nx,ny,nz,shb,istatus)
      call get_modelfg_3d(i4time_sys_drop,'OM ',nx,ny,nz,omb,istatus)
      if(istatus.ne.1)then
         print*,'Error returned: get_modelfg_3d'
      endif
      call get_modelfg_2d(i4time_sys_drop,'PSF',nx,ny,psb,istatus)
      if(istatus.ne.1)then
         print*,'Error returned: get_modelfg_2d'
      endif
c
c use this background with the input background to advance the analysis
c to time i4time_sys_drop.
c
      lapsu=lapsu+(ub-ub_a)
      lapsv=lapsv+(vb-vb_a)
      lapstemp=lapstemp+(tb-tb_a)
      lapssh=lapssh+(shb-shb_a)
      lapsphi=lapsphi+(phib-phib_a)
      lapsomo=lapsomo+(omb-omb_a)
      ps=ps+(psb-psb_a)

c return backgrounds at i4time_sys_drop.
      ub_a=ub
      vb_a=vb
      tb_a=tb
      shb_a=shb
      phib_a=phib
      omb_a=omb
      psb_a=psb

      deallocate (ub,vb,tb,phib,shb,omb,psb)
      
      return
      end
c
c-------------------------------------------------------
c
      subroutine read_wind3d_wgi(rstats,istatus)

      implicit  none

      character a9_time*9
      character directory*255
      character filename*255
      character logdir*255
      character dum1*255
      character var(20)*20
      integer   cnt,vnum,vnumi
      integer   i,is,ie,len,istatus,i4time_sys
      logical   found_line,lexist
      real      rstats(7)

      call get_directory('log',logdir,len)
      call get_systime(i4time_sys,a9_time,istatus)
      if(istatus.ne.1)return

      filename=logdir(1:len)//'/wind.wgi.'//a9_time
      call s_len(filename,len)
      print*,'filename = ',filename(1:len)
      inquire(file=filename,exist=lexist)
      found_line=.false.
      if(lexist)then
        open (11, file=filename,form='formatted',status='old',err=50) 
        Do while (.not.found_line)
         read(11,100,end=1)dum1
         if(dum1(1:21).eq.'Obs minus First Guess')then
            read(11,*,err=1)  !this one to read the column labels
            do while (.not.found_line)
               read(11,100,err=1)dum1
               if(dum1(1:6).eq."   PIN".or.dum1(1:8).eq."  DROPSN")then
                  found_line=.true.
               endif
            enddo
         endif
        enddo
      else
        print*,'file not found: ',filename(1:len)
        istatus = 0
        return
      endif
100   format(a)

      i=255
      do while (i.gt.0)
         i=i-1
         if(dum1(i:i).ne.' ')then
            ie=i
            i=0
         endif
      enddo

      var=' '
      i=7
      if(dum1(1:8).eq."  DROPSN")i=9
      do while (i.ne.0)
         if(dum1(i:i).ne.' ')then
            is=i
            i=0
         else
            i=i+1
         endif
      enddo

      dum1=dum1(is:ie)
      ie=ie-is+1
      vnumi=0
      vnum=1
      cnt=0
      do i=1,ie
         if(dum1(i:i) .ne. ' ')then
            cnt=cnt+1
            if(vnum.ne.vnumi)vnumi=vnumi+1
            var(vnumi)(cnt:cnt)=dum1(i:i)
         else
            cnt=0
            if(vnum.eq.vnumi)vnum=vnum+1
         endif
      enddo

c     do i=1,vnumi
c        print*,'wgi num ',i,' = ',var(i)
c     enddo

      do i=1,vnumi
         read(var(i),'(f6.2)')rstats(i)
      enddo

c     do i=1,vnumi
c        print*,'wgi num ',i,' = ',var(i)
c     enddo

      return

50    print*, 'Error opening file: ',filename
1     print*, 'Suspect read in the wgi file'
      istatus=0
      return
      end
