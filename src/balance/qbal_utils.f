
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
      terscl=sqrt(sumt/(nx*ny))

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

      subroutine advance_grids(i4time_sys,nx,ny,nz,ub,vb,tb,phib
     .,shb,omb,lapsphi,lapstemp,lapsu,lapsv,lapssh,omo,ps,istatus)

      implicit none

      include 'bgdata.inc'

      integer nx,ny,nz

      real    phib(nx,ny,nz)
     .       ,tb(nx,ny,nz)
     .       ,ub(nx,ny,nz)
     .       ,vb(nx,ny,nz)
     .       ,shb(nx,ny,nz)
     .       ,omb(nx,ny,nz)

      real    lapsu(nx,ny,nz)
     .       ,lapsv(nx,ny,nz)
     .       ,lapssh(nx,ny,nz)
     .       ,lapstemp(nx,ny,nz)
     .       ,lapsphi(nx,ny,nz)
     .       ,omo(nx,ny,nz)
     .       ,ps(nx,ny)

      character*256 bgpaths(maxbgmodels)
      character*132 cmodels(maxbgmodels)
      integer   bgmodels(maxbgmodels)
      integer   itime_inc

      character(len=256), allocatable :: names(:)
      character(len=256), allocatable :: reject_names(:)
      character(len=256), allocatable :: bg_names(:)

      integer   max_files
      parameter (max_files=200)
      integer accepted_files
      integer bg_files
      integer oldest_forecast
      integer forecast_length
      integer max_forecast_delta
      integer rejected_cnt
      logical use_analysis

      integer i4time_sys
     .       ,lenm
     .       ,istatus

      print*,'Test subroutine advance_grid: doesnt do much yet.'

c     1. use readbgdata with input (advanced) i4time_sys and
c        path to do bgdata to read LAPS_FUA. Because this routine
c        expects a full path/filename, we must determine this apriori
c        which is a lga/get_acceptable_files function.

      call get_background_info(bgpaths,bgmodels,oldest_forecast
     +,max_forecast_delta,forecast_length,use_analysis,cmodels
     +,itime_inc)

      call s_len(cmodels(1),lenm)

      if(cmodels(1).ne.'LAPS_FUA')then
         print*,' Error: analysis and background grids are'
     .,' advanced in time only when background model is LAPS_FUA'
     .,' as specified in background.nl'
            print*,' Current background: ',cmodels(1)(1:lenm)
         return
      endif

      allocate (names(max_files),reject_names(max_files)
     .         ,bg_names(max_files))

!     call get_acceptable_files(i4time_sys,bgpaths(1),bgmodels(1)
!    +        ,names,max_files,oldest_forecast,max_forecast_delta
!    +        ,use_analysis,bg_files,accepted_files,forecast_length
!    +        ,cmodels(1),nx,ny,nz,reject_names,rejected_cnt)

      deallocate (names,reject_names,bg_names)

      return
      end
