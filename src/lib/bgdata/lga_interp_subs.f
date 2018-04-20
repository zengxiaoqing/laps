cdis    Forecast Systems Laboratory
cdis    NOAA/OAR/ERL/FSL
cdis    325 Broadway
cdis    Boulder, CO     80303
cdis 
cdis    Forecast Research Division
cdis    Local Analysis and Prediction Branch
cdis    LAPS 
cdis 
cdis    This software and its documentation are in the public domain and 
cdis    are furnished "as is."  The United States government, its 
cdis    instrumentalities, officers, employees, and agents make no 
cdis    warranty, express or implied, as to the usefulness of the software 
cdis    and documentation for any purpose.  They assume no responsibility 
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis    
cdis    Permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  All modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making 
cdis    the modifications.  If significant modifications or enhancements 
cdis    are made to this software, the FSL Software Policy Manager  
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
c
c===============================================================================
c
      subroutine time_interp(dir,ext,nx,ny,nz,ngrids,pr,ht_1d,
     .                  i4time_valid1,i4time_valid2,i4time_now,
     .                  time1,fcst1,time2,fcst2)

      use mem_namelist, ONLY: vertical_grid
c
      implicit none
      include 'bgdata.inc'
c
      integer nx,ny,nz
      integer n,ngrids,ngrids_loc
      integer warncnt
c
      integer   i4time_valid1,
     .          i4time_valid2,
     .          i4time_now,
     .          i4time_needed_in,
     .          i4time_nearest,
     .          time1,time2,
     .          fcst1,fcst2,
     .          ip(nz),lend,
     .          newfcst,
     .          imin,ihour,
     .          i,j,k,kk,
     .          istatus,nstatus,ishow_timer

      logical   lexist
c
      real   pr(nz),ht_1d(nz),weight

      real,  allocatable :: 
     .       grid1(:,:,:),
     .       grid2(:,:,:),
     .       gridn(:,:,:)
c
      integer        nan_flag

      character*13   cfilename
      character*255  cfilespec
      character*(*)  dir
      character*(*)  ext
      character*3    var(nz,ngrids)
      character*4    lvl_coord(nz)
      character*10   units(nz)
      character*125  comment(nz)
      character*9    fname9
      character*4    af

!!    add by Huiling Yuan  for getting background information, 20120905, AAA001
      character*256 bgpaths(maxbgmodels)
      character*132 cmodel(maxbgmodels)
      integer bgmodels(maxbgmodels)
      integer forecast_length
      integer itime_inc,ntmin,ntmax
      integer nbgm
      logical sfc_bkgd
      logical use_analysis, use_forecast
      logical smooth_fields
      logical lgb_only
      integer  i4reftime,i4valtime
!!    end by Huiling Yuan, AAA001
c_______________________________________________________________________________
c
      write(6,*)' Subroutine time_interp...',trim(ext),ngrids

!! add Huiling Yuan 20120905, AAA002
c
c Read information from static/background.nl,  from src/lib/laps_io.f 
c
      call get_background_info(bgpaths,bgmodels
     +,forecast_length
     +,use_analysis,use_forecast
     +,cmodel,itime_inc,smooth_fields,sfc_bkgd,ntmin,ntmax,lgb_only)

      print*, '----------------------------------------------'
      print*, 'time_interp: use_analysis = ',use_analysis
      print*, '             use_forecast = ',use_forecast
      print*, '----------------------------------------------'
!! end Huiling Yuan AAA002

      if(ext.eq.'lga') then
         do k=1,nz
            if(vertical_grid .ne. 'SIGMA_HT')then
               ip(k)=int(pr(k))
               var(k,1)='HT '
            else
               ip(k)= int(ht_1d(k)) 
               var(k,1)='P3 '
            endif
            var(k,2)='T3 '
            var(k,3)='SH '
            var(k,4)='U3 '
            var(k,5)='V3 '
            if(vertical_grid .ne. 'SIGMA_HT')then
              var(k,6)='OM '
            else
              var(k,6)='W3 '
            endif
         enddo
      else
         do k=1,nz
            ip(k)=0
         enddo
         var(1,1)='USF'
         var(1,2)='VSF'
         var(1,3)='TSF'
         var(1,4)='PSF'
         var(1,5)='SLP'
         var(1,6)='RSF'
         var(1,7)='DSF'
         var(1,8)='P  '
         var(1,9)='TGD'
         var(1,10)='R01'
      endif


      weight=float(i4time_valid2-i4time_now)/
     .       float(i4time_valid2-i4time_valid1)
      print*,'Time interp weight = ',weight

      write(af,'(i4.4)') fcst1/3600
      call  make_fnam_lp(time1, fname9, nstatus)
      print*,'Time interp file 1: ',fname9,af,'.'//ext
      write(af,'(i4.4)') fcst2/3600
      call  make_fnam_lp(time2, fname9, nstatus)
      print*,'Time interp file 2: ',fname9,af,'.'//ext

!     add Huiling Yuan, 20120906, AAA003
      if((use_analysis. eqv. .true.) .and. 
     1   (use_forecast .eqv. .false.))then   
      newfcst=i4time_now-i4time_valid1
      call make_fnam_lp(time2,fname9,istatus)
      imin=mod(newfcst,3600)/60
      ihour=newfcst/3600
      write(af,'(2i2.2)') ihour,imin
      print *,'Time interp output file: ',fname9//af,'.'//ext(1:3)

      else
      newfcst=fcst1-(i4time_valid2-i4time_now)
      call make_fnam_lp(time1,fname9,istatus)
      imin=mod(newfcst,3600)/60
      ihour=newfcst/3600
      write(af,'(2i2.2)') ihour,imin
      print *,'Time interp output file: ',fname9//af,'.'//ext(1:3)
      endif   ! end if block use_analysis=true, Huiling Yuan, 20120906, AAA003

      call s_len(dir,lend)
      cfilespec=dir(1:lend)//'/'//fname9//af//'.'//ext
      inquire(file=cfilespec,exist=lexist)
      if(lexist)then
         print*,'Output file already exists in ',ext
         print*,'Returning to lga_driver. No time interp.'
         return
      endif

      allocate (grid1(nx,ny,nz)
     &         ,grid2(nx,ny,nz)
     &         ,gridn(nx,ny,nz))


      do n=1,ngrids

         call read_laps(time1,time1+fcst1,dir,ext,
     .        nx,ny,nz,nz,var(1,n),ip,lvl_coord,units,comment,
     .        grid1,istatus)
c
         if(istatus.ne.1) then
            print *, 'ERROR returned from read_laps, time1'
            print *, 'n / ngrids = ',n,ngrids                        
            stop 'time_interp'
         endif

         call read_laps(time2,time2+fcst2,dir,ext,
     .        nx,ny,nz,nz,var(1,n),ip,lvl_coord,units,comment,
     .        grid2,istatus)
         if(istatus.ne.1) then
            print *, 'ERROR returned from read_laps, time2'
            print *, 'n / ngrids = ',n,ngrids                        
            stop 'time_interp'
         endif
c
c *** Do interpolation with time for each new file.
c
         warncnt = 0

         call check_nan3(grid1,nx,ny,nz,nan_flag)
         if(nan_flag .ne. 1) then
            print *,' ERROR: NaN found in grid1 array '
            goto 99
         endif
c
         call check_nan3(grid2,nx,ny,nz,nan_flag)
         if(nan_flag .ne. 1) then
            print *,' ERROR: NaN found in grid2 array '
            goto 99
         endif

         do k=1,nz
         do j=1,ny
         do i=1,nx
            gridn(i,j,k)= (1.- weight)*grid1(i,j,k) +
     +                               weight *grid2(i,j,k)

         enddo
         enddo
         enddo

!        add Huiling Yuan, 20120906, AAA004
         if((use_analysis. eqv. .true.) .AND. 
     1      (use_forecast .eqv. .false.)      )then   
            i4reftime=time2
            i4valtime=time2+newfcst
         else
            i4reftime=time1
            i4valtime=time1+newfcst
         endif    ! end Huiling Yuan, AAA004
 
c
!        istatus = ishow_timer()
         comment(1) = 'Time Interpolated: '//comment(1)
         if(ext.eq.'lga')then
!YHL 20120906            call write_laps(time1,time1+newfcst,dir,ext,
            call write_laps(i4reftime,i4valtime,dir,ext,
     .           nx,ny,nz,nz,var(1,n),ip,lvl_coord,units,comment,
     .           gridn,istatus)
         else
!YHL 20120906            call write_laps(time1,time1+newfcst,dir,ext,
            call write_laps(i4reftime,i4valtime,dir,ext,
     .           nx,ny,1,1,var(1,n),ip,lvl_coord,units,comment,
     .           gridn,istatus)
         endif

      enddo

99    deallocate (grid1, grid2, gridn)
c
      return
      end

c     subroutine erase_file(inittime,validtime,dir,ext)
c     integer inittime,validtime, istatus, rename
c     character*(*) dir, ext
c     character*256 filename
c     character*13 fname
c     call make_fnam13_lp(inittime,validtime, fname, istatus)
c     print*,inittime,validtime, fname
c     call s_len(dir,len_dir)
c     if(dir(len_dir:len_dir) .ne. '/') then
c        dir(len_dir:len_dir)='/'
c        len_dir=len_dir+1
c     endif
c     write(filename,*) dir(1:len_dir)//fname//'.'
c     istatus = rename(filename(1:len_dir+15)//ext(1:3)
c    +     ,filename(1:len_dir+15)//'bad')
c     print*,'renamed ',filename(1:len_dir+15)//ext(1:3),
c    +     ' with istatus= ',istatus
c     return
c     end
