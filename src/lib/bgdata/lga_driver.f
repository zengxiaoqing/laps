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
      subroutine lga_driver(nx_laps,ny_laps,nz_laps,luse_sfc_bkgd
     .    ,mode
     .    ,laps_cycle_time,bgmodel,bgpath,cmodel,reject_cnt
     .    ,reject_names,bg_names,max_files,accepted_files
     .    ,n_written,c_ftimes_written
     .    ,i4time_now,smooth_fields,lgb_only,lga_status)

c KML: CHANGES MADE APRIL 2004
c tdbg_sfc (model 2m dew point) is now being read in from subroutine read_bgdata
c tdbg_sfc is then horizontally interpolated to LAPS grid (td_sfc)
c ht_sfc and td_sfc are passed into subroutine sfcbkgd -> now called sfcbkgd_sfc (JRS)
c KML: END

 
      use mem_namelist, ONLY: vertical_grid
      use storage_module, ONLY: get_plvls
 
      implicit none
      integer, parameter :: maxbglvl = 52 

      real    prtop ! pa
      include 'bgdata.inc'
      real     badflag
      include 'laps_sfc.inc'
c
      integer nx_laps,ny_laps,nz_laps,     !LAPS grid dimensions
     .        nx_pr  ,ny_pr,               !LAPS pressure grid dimensions
     .        nx_bg  ,ny_bg,               !Background model grid dimensions
     .        nzbg_ht,nzbg_tp,nzbg_sh,
     .        nzbg_uv,nzbg_ww,
     .        nlvl_grib,                   !Actual number of grib levels
     .        max_files, lga_status,
     .        laps_cycle_time,
     .        bgmodel,nbg,lncf,mode

      integer   ishow_timer
      integer   init_timer
      integer   itstatus(10)
      integer   itstatus_rot,istat_alloc
      integer   icnt
      integer   igrx,igry
      integer ixmin, ixmax, iymin, iymax, i_perim
      integer i_mx, i_mn, j_mx, j_mn, nan_flag
      real diff, diff_mx, diff_mn
      real xmin, xmax, ymin, ymax, bgres, grid_spacing_cen_m

      real Lon0,Lat0,Lat1,dlat,dlon
      real sw(2),ne(2)
      real cenlat,cenlon
      real dx,dy
      real psatoz

      character*256 bgpath
      character*256 bg_names(max_files)         ! should be basenames
      character*256 reject_names(max_files)
      character*132 cmodel
      character*150 static_dir,filename

      integer warncnt
      
      integer n_written, iw
      character*14 c_ftimes_written(100), a14_time

c sfc namelist stuff. for reduced pressure calc
      integer use_lso_qc, skip_internal_qc, itheta
      logical l_require_lso, luse_sfc_bkgd, ltest_vertical_grid
      logical lgb_only
      real    redp_lvl,del,gam,ak
      real    bad_t,bad_td,bad_u,bad_v,bad_p
      real    bad_mp,bad_th,bad_the
      real    bad_vis,bad_tb8
      real    thresh_t,thresh_td,thresh_mslp
      real    rms_wind
c
c *** Background model grid data.
c
c
c *** sfc background arrays.
c
      real, allocatable  :: prbg_sfc(:,:)    !In Pascals
      real, allocatable  :: uwbg_sfc(:,:)
      real, allocatable  :: vwbg_sfc(:,:)
      real, allocatable  :: shbg_sfc(:,:)
      real, allocatable  :: tdbg_sfc(:,:)
      real, allocatable  :: tpbg_sfc(:,:)
      real, allocatable  :: t_at_sfc(:,:)
      real, allocatable  :: htbg_sfc(:,:)
      real, allocatable  :: mslpbg(:,:)
      real, allocatable  :: pcpbg(:,:)       !Precip at surface, ACPC (k/m^2)
      real, allocatable  :: crefbg(:,:)      !Composite Reflectivity (dBZ)
      real, allocatable  :: tpwbg(:,:)       !Total Precipitable Water
      real, allocatable  :: cwatbg(:,:)      !Cloud Water (Vert Integrated)
      real, allocatable  :: swibg(:,:)       !Downward Solar Radiation (GHI)

c
c *** 3D background arrays.
c
      real, allocatable  :: prbght(:,:,:)    !Pressure (mb) height levels
      real, allocatable  :: prbgsh(:,:,:)    !Pressure (mb) humidity levels
      real, allocatable  :: prbguv(:,:,:)    !Pressure (mb) u/v components
      real, allocatable  :: prbgww(:,:,:)    !Pressure (mb) omega levels
      real, allocatable  :: htbg(:,:,:)      !Background heights (m)
      real, allocatable  :: tpbg(:,:,:)      !Background temps (K)
      real, allocatable  :: shbg(:,:,:)      !Background humidity (gm/m3)
      real, allocatable  :: uwbg(:,:,:)      !Background u component (m/s)
      real, allocatable  :: vwbg(:,:,:)      !Background v component (m/s)
      real, allocatable  :: wwbg(:,:,:)      !Background omega (m/s)
      real plvl_grib(maxbglvl) ! Dimension is maxbglvl in 'degrib_nav' routine
c
c *** Intermediate arrays for background data vertically
c     interpolated to LAPS isobaric levels (on the model horizontal grid).
c
      real, allocatable :: htvi(:,:,:)
      real, allocatable :: tpvi(:,:,:)
      real, allocatable :: shvi(:,:,:)
      real, allocatable :: uwvi(:,:,:)
      real, allocatable :: vwvi(:,:,:)
      real, allocatable :: wwvi(:,:,:)
      real, allocatable :: prvi(:,:,:) ! Pressure (mb)                 

      integer, allocatable :: msgpt(:,:)
c
c *** Background data interpolated to LAPS grid.
c
      real      ht(nx_laps,ny_laps,nz_laps), !Height (m)
     .          tp(nx_laps,ny_laps,nz_laps), !Temperature (K)
     .          sh(nx_laps,ny_laps,nz_laps), !Specific humidity (kg/kg)
     .          uw(nx_laps,ny_laps,nz_laps), !!U-wind (m/s)
     .          vw(nx_laps,ny_laps,nz_laps), !V-wind (m/s)
     .          ww(nx_laps,ny_laps,nz_laps), !W-wind (pa/s)
     .          pr1d_pa(nz_laps),            !LAPS pressures (pa)
     .          pr1d_mb(nz_laps),            !LAPS pressures (mb)
     .          sigma1d(nz_laps),            !LAPS vert grid (SIGMA_P grid only)
     .          ht_1d(nz_laps),              !LAPS vert grid (SIGMA_HT grid only)
     .          lat(nx_laps,ny_laps),        !LAPS lat
     .          lon(nx_laps,ny_laps),        !LAPS lon
     .          topo(nx_laps,ny_laps),       !LAPS avg terrain
     .          grx(nx_laps,ny_laps),        !hinterp factor
     .          gry(nx_laps,ny_laps),        !hinterp factor
     .          ht_sfc(nx_laps,ny_laps),     !first guess terrain
     .          td_sfc(nx_laps,ny_laps),     !2m dewpoint
     .          td_sfc_hi(nx_laps,ny_laps),  !2m dewpoint (on hi-res sfc)
     .          tp_sfc(nx_laps,ny_laps),     !2m surface temperature
     .          t_sfc (nx_laps,ny_laps),     !sea/land surface temp
     .          sh_sfc(nx_laps,ny_laps),
     .          qsfc(nx_laps,ny_laps),
     .          uw_sfc(nx_laps,ny_laps),
     .          vw_sfc(nx_laps,ny_laps),
     .          pr_sfc(nx_laps,ny_laps),     !Stn pressure
     .          rp_sfc(nx_laps,ny_laps),     !Reduced pressure
     .          mslp(nx_laps,ny_laps),       !in Pascals
     .          pcp_sfc(nx_laps,ny_laps),
     .          cw_sfc(nx_laps,ny_laps),
     .          dum1_2d(nx_laps,ny_laps),
     .          dum2_2d(nx_laps,ny_laps)

      real, allocatable :: rp_lvl(:,:)       !Reduced pressure lvl
      real, allocatable :: rp_tp(:,:)        !Reduced pressure temp (holder)
      real, allocatable :: rp_sh(:,:)        !Reduced pressure sh   (holder)
      real, allocatable :: rp_td(:,:)        !Reduced pressure td   (holder)
      real, allocatable :: prgd_pa(:,:,:)    !Pressure 3D/1D on LAPS Grid (Pa)
c
      real      ssh2,                        !Function name
     .          shsat,cti,
     .          htave,tpave,shave,uwave,vwave,
     .          rmaxvv,rminvv,
     .          tp_sfc_c,td_sfc_c
c
      integer   ct,  reject_cnt,
     .          ihour,imin,
     .          i4time_now,
     .          lga_files,lga_times(max_files),
     .          lga_valid,
     .          bg_times(max_files),accepted_files, bglen,
     .          i4time_bg_valid(max_files),
     .          bg_valid(max_files),
     .          valid_bg(max_files),time_bg(max_files),
     .          bgvalid,
     .          i,ic,ii,j,jj,k,kk,l,ldl,lf,
     .          istatus

      integer   istatus_prep(max_files)

      integer   i4lgatime(max_files)
      integer   i4bgtime
      integer   nlga
c
      character*256 lgapath
      character*256 lga_names(max_files)
      character*256 names(max_files)
      character*256 fname_bg(max_files)

      character*13  fname13,fname9_to_wfo_fname13
      character*13  cvt_fname13_to_wfo_fname13
      character*13  cvt_wfo_fname13_to_fname13
      character*9   wfo_fname13_to_fname9

      character*6   c6_maproj
      character*2   gproj
      character*1   cgrddef
      character*200 fullname
      character*256 outdir
      character*31  ext
      character*3   c_ext
      character*125 comment(nz_laps)
      character*4   af_bg(max_files)
      character*10  c_domain_name
      character*200 c_dataroot
      character*256 cfname

      integer len_dir,ntime, nf
c
      logical llapsfua
      logical linterp
      logical smooth_fields
      logical lgrid_missing
      logical lhif_tsfc
      logical wrapped, l_bilinear

      data ntime/0/
      data ext/'lga'/

      interface
c
       subroutine read_bgdata(nx_bg,ny_bg
     +,nzbg_ht,nzbg_tp,nzbg_sh,nzbg_uv,nzbg_ww,ctype
     +,bgpath,fname_bg,af_bg,fullname,cmodel,bgmodel
     +,prbght,prbgsh,prbguv,prbgww
     +,htbg,tpbg,uwbg,vwbg,shbg,wwbg
     +,htbg_sfc,prbg_sfc,shbg_sfc,tdbg_sfc,tpbg_sfc
     +,t_at_sfc,uwbg_sfc,vwbg_sfc,mslpbg,pcpbg,crefbg,tpwbg,cwatbg,swibg
     +,istatus)
c
         real  :: prbg_sfc(nx_bg,ny_bg)
         real  :: uwbg_sfc(nx_bg,ny_bg)
         real  :: vwbg_sfc(nx_bg,ny_bg)
         real  :: shbg_sfc(nx_bg,ny_bg)
         real  :: tdbg_sfc(nx_bg,ny_bg)
         real  :: tpbg_sfc(nx_bg,ny_bg)
         real  :: t_at_sfc(nx_bg,ny_bg)
         real  :: htbg_sfc(nx_bg,ny_bg)
         real  :: mslpbg(nx_bg,ny_bg)
         real  :: pcpbg(nx_bg,ny_bg)
         real  :: crefbg(nx_bg,ny_bg)
         real  :: tpwbg(nx_bg,ny_bg)
         real  :: cwatbg(nx_bg,ny_bg)
         real  :: swibg(nx_bg,ny_bg)
c
         real  :: prbght(nx_bg,ny_bg,nzbg_ht)
         real  :: prbgsh(nx_bg,ny_bg,nzbg_sh)
         real  :: prbguv(nx_bg,ny_bg,nzbg_uv)
         real  :: prbgww(nx_bg,ny_bg,nzbg_ww)
         real  :: htbg(nx_bg,ny_bg,nzbg_ht)
         real  :: tpbg(nx_bg,ny_bg,nzbg_tp)
         real  :: shbg(nx_bg,ny_bg,nzbg_sh)
         real  :: uwbg(nx_bg,ny_bg,nzbg_uv)
         real  :: vwbg(nx_bg,ny_bg,nzbg_uv)
         real  :: wwbg(nx_bg,ny_bg,nzbg_ww)

         character*4   af_bg
         character*5   ctype
         character*132 cmodel
         character*256 bgpath
         character*256 fname_bg
         character*200 fullname
         integer       bgmodel
         integer       nx_bg
         integer       ny_bg
         integer       nzbg_ht
         integer       nzbg_tp
         integer       nzbg_sh
         integer       nzbg_uv
         integer       nzbg_ww
         integer       istatus
       end subroutine

       subroutine get_bkgd_mdl_info(bgmodel
     &,cmodel,fullname,nx,ny,nzbg_ht,nzbg_tp
     &,nzbg_sh,nzbg_uv,nzbg_ww
     &,gproj,dlat,dlon,centrallat,centrallon,dx,dy
     &,Lat0,Lat1,Lon0,sw,ne,cgrddef,istatus)

        character*200 fullname
        character*132 cmodel
        character*2   gproj
        character*1   cgrddef
        integer       istatus
        integer       bgmodel
        integer       nx,ny
        integer       nzbg_ht
        integer       nzbg_tp
        integer       nzbg_sh
        integer       nzbg_uv
        integer       nzbg_ww
        real          Lat1
        real          Lat0,Lon0
        real          sw(2),ne(2)
        real          centrallat
        real          centrallon
        real          dx,dy

       end subroutine

       subroutine init_hinterp(nx_bg,ny_bg,nx_laps,ny_laps,gproj,
     .     lat,lon,grx,gry,bgmodel,cmodel,wrapped)

         character  gproj*2
         character  cmodel*132
         integer nx_bg,ny_bg,nx_laps,ny_laps,bgmodel
         logical wrapped
         real   lat(nx_laps,ny_laps)
     .         ,lon(nx_laps,ny_laps)
     .         ,grx(nx_laps,ny_laps)
     .         ,gry(nx_laps,ny_laps)

       end subroutine

       subroutine filter_2dx(field,ix,iy,iz,smth)
         integer ix,iy,iz
         real field(ix,iy,iz)
         real smth
       end subroutine


      end interface

      warncnt=0 
c_______________________________________________________________________________
c *** Get LAPS lat, lons.
c
      print *,'in lga_driver'
      linterp=.true.
      lga_status=0

      call s_len(cmodel,ic)

      call s_len(bgpath,bglen)
      if(bgpath(bglen:bglen).ne.'/')then
         bglen=bglen+1
         bgpath(bglen:bglen)='/'
      endif
      cfname=bgpath(1:bglen)//bg_names(1)
      call s_len(cfname,lncf)

      write(6,*)' bgpath = ',bgpath(1:bglen)
      write(6,*)' bg_names(1) = ',TRIM(bg_names(1))
      write(6,*)' lncf = ',lncf
      write(6,*)' cfname(1:lncf) = ',cfname(1:lncf)

C WNI Add a section to identify wrapped grid and set the wrapped flag
      wrapped = .FALSE.              ! WNI
      IF (   bgmodel .eq. 6 .or.     ! WNI
     .       bgmodel .eq. 8 .or.     ! WNI
     .      (bgmodel .eq. 13 .and. cmodel(1:ic) .eq. 'GFS')   .or.  ! SCA
     .      (bgmodel .eq. 13 .and. cmodel(1:ic) .eq. 'ECMWF') .or.  ! SCA
     .      (bgmodel .eq. 13 .and. cmodel(1:ic) .eq. 'FIM')   .or.  ! SCA
     .       bgmodel .eq. 10.) THEN  ! WNI
        wrapped = .true.             ! WNI
      ENDIF                          ! WNI
C WNI END ADDITON

      if(lgb_only .and. (.not. wrapped))then
          l_bilinear = .true.
      else
          l_bilinear = .false.
      endif

      istatus=ishow_timer()

      print *,'calling get_bkgd_mdl_info: cfname = ',cfname

      call get_bkgd_mdl_info(bgmodel,cmodel,cfname
     &,nx_bg,ny_bg,nzbg_ht,nzbg_tp,nzbg_sh,nzbg_uv
     &,nzbg_ww ,gproj,dlat,dlon,cenlat,cenlon,dx,dy
     &,Lat0,Lat1,Lon0,sw,ne,cgrddef,istatus)

      if(istatus.ne.1)then
         print*,'Error getting background model information'
         ! Yuanfu: Check istatus and show hint for users
         if (istatus .eq. -4) 
     &     print*,'Check background.nl for bgmodel and cmodel'
         print*,'Returning from lga_driver'
         return
      endif

      istatus=ishow_timer()
      write(6,*)' Returned from get_bkgd_mdl_info'

      llapsfua=.false.
      if(bgmodel.eq.0)then
         if(cmodel(1:ic).eq.'MODEL_FUA'.or.
     &      cmodel(1:ic).eq.'LAPS_FUA' .or.
     &      cmodel(1:ic).eq.'LAPS')then
            llapsfua=.true.
            if(cmodel(1:ic).eq.'LAPS_FUA'.or.
     &         cmodel(1:ic).eq.'LAPS')then

               linterp = .false.

               print*,'*************************************'
               print*,'Interpolation to domain not necessary'
               print*,'*************************************'

            endif
         endif
      endif

      call init_gridconv_cmn(gproj,nx_bg,ny_bg,nzbg_ht
     &,dlat,dlon,cenlat,cenlon,Lat0,Lat1,Lon0
     &,sw(1),sw(2),ne(1),ne(2),cgrddef,istatus)

      print *
      print *, ' Background information: '
      print *, '-------------------------'
      print *, 'nx/ny:              ',nx_bg,ny_bg
      print *, 'nz(ht,tp,sh,uv,ww): '
     &,nzbg_ht,nzbg_tp,nzbg_sh,nzbg_uv,nzbg_ww
      print *, 'bgpath:             ',bgpath(1:bglen)
      print *, 'bgmodel:            ',bgmodel
      print *, 'cmodel:             ',cmodel(1:ic)
      print *

      call get_directory('static',outdir,len_dir)    

      call find_domain_name(c_dataroot,c_domain_name,istatus)
      if(istatus.ne.1)then
          print*,'Error returned: find_domain_name'
          return
      endif

      call get_laps_domain(nx_laps,ny_laps,c_domain_name
     +,lat,lon,topo,istatus)
      if (istatus.lt.1)then
          print *,'Error reading lat, lon, topo data.'
          stop
      endif
c
c *** Specify model path, extension for write laps routine.
c
      call get_directory('lga',outdir,len_dir)
      print *,'output dir will be ',outdir(1:len_dir)
c
c *** get LAPS pressure OR height levels.  Using pressures.nl / heights.nl
c
      write(6,*)' Vertical grid is: ',vertical_grid
      if(vertical_grid .eq. 'PRESSURE')then ! PRESSURE
          print*,'get 1d pressures'
          call get_pres_1d(i4time_now,nz_laps,pr1d_pa,istatus)
          if(istatus.ne.1)then
             print*,'Error returned from get_pres_1d'
             print*,'Check pressures.nl or nk_laps in nest7grid.parms'
             stop
          endif
          pr1d_mb(:)=pr1d_pa(:)/100.  ! Pa to mb

          nx_pr = 1
          ny_pr = 1
!         allocate (prgd_pa(nx_pr,ny_pr,nz_laps), STAT=istat_alloc)
!         if(istat_alloc .ne. 0)then
!             write(6,*)' ERROR: Could not allocate prgd_pa'
!             stop
!         endif
!         prgd_pa(1,1,:)=pr1d_pa(:)   ! 1D Pa

      elseif(vertical_grid .eq. 'SIGMA_HT')then
          nx_pr = nx_laps
          ny_pr = ny_laps
          allocate (prgd_pa(nx_pr,ny_pr,nz_laps))

      elseif(vertical_grid .eq. 'SIGMA_P')then
          if(.not. luse_sfc_bkgd)then ! test for valid setting
             write(6,*)' ERROR: luse_sfc_bkgd is FALSE in namelist'
             write(6,*)' This should be TRUE for a SIGMA_P grid'
             stop
          endif

          print*,'get 1d sigmas'
          call get_sigma_1d(nz_laps,sigma1d,istatus)
          if(istatus.ne.1)then
             print*,'Error returned from get_sigma_1d'
             print*,'Check sigmas.nl or nk_laps in nest7grid.parms'
             stop
          endif
          nx_pr = nx_laps
          ny_pr = ny_laps

          allocate (prgd_pa(nx_laps,ny_laps,nz_laps))

      endif
c
c *** Determine which of the "bg_names" has not already been processed
c
      do j=1,max_files
         names(j)=bg_names(j)
      enddo

      do j=1,accepted_files
         call i4time_fname_lp(names(j)(1:9),bg_times(j),istatus)
         imin=0
         if(llapsfua)then ! may need > 100 hour support
            read(bg_names(j)(10:11),'(i2)')ihour
            read(bg_names(j)(12:13),'(i2)')imin
         else
            if(cmodel(1:ic) .eq. 'HRRR')then    
               read(bg_names(j)(10:11),'(i2)')ihour
               read(bg_names(j)(12:13),'(i2)')imin       
            elseif(cmodel(1:ic) .eq. 'UM' .OR.
     1             cmodel(1:ic) .eq. 'GUM'     )then
               read(bg_names(j)(10:11),'(i2)')ihour
            else
               read(bg_names(j)(12:13),'(i2)')ihour
            endif
         endif
         bg_valid(j)=ihour*3600+imin*60
         i4time_bg_valid(j)=bg_times(j)+bg_valid(j)
      enddo

      c_ext='lga'
      if(lgb_only)c_ext='lgb'
      call get_directory(c_ext,lgapath,ldl)
      call get_file_times(lgapath,max_files,lga_names
     1                      ,lga_times,nlga,istatus)

      if(nlga.gt.0)then
         k=1
         do while (k.le.nlga)

            read(lga_names(k)(ldl+10:ldl+11),'(i2)')ihour
            read(lga_names(k)(ldl+12:ldl+13),'(i2)')imin

            lga_valid=ihour*3600+imin*60

            do j=1,accepted_files

               if(bg_times(j).eq.lga_times(k).and.
     .            bg_valid(j).eq.lga_valid)      then
c
c ok, since we've processed it, lets add it to the reject list
c
                  reject_cnt=reject_cnt+1
                  reject_names(reject_cnt)=names(j)
c                 names(j)=' '           ! commented out 12/2009 by SA

                  print*,'reject_cnt/reject_names'
                  print*,'cnt/time: '
     +                  ,reject_names(reject_cnt)

                  print *,'NOTE that LGA/LGB file exists: ',lga_names(k)
!                 lga_status = 1
!                 goto 900

               endif
            enddo
            k=k+1
         enddo
      endif

      nbg=0
      do 33 j=1,accepted_files

         if(names(j).ne.' ')then

          do k=1,reject_cnt
           if(reject_names(k).eq.names(j))goto 33 !then
c             names(j)=' '
c          endif

          enddo

          nbg = nbg+1
          if(bgmodel.eq.4.or.bgmodel.eq.10)then
             fname_bg(nbg)=
     1 cvt_fname13_to_wfo_fname13(bg_names(j)(1:13))
          else
             fname_bg(nbg)=bg_names(j)
          endif

          af_bg(nbg)=bg_names(j)(10:13)
          time_bg(nbg)=bg_times(j)
          valid_bg(nbg)=bg_valid(j)

         endif

33    continue

      if(nbg.eq.0)then
         print*,'No new model background to process'
         lga_status = 1
      endif
c
c ****** Read background model data.
c
c     print*,'process new model background'

      do nf=1,nbg
 
       bgvalid=time_bg(nf)+valid_bg(nf)

!      Check a14_time and compare to previously processed files
       call make_fnam_lp(time_bg(nf),a14_time(1:9),istatus)       
       call make_fcst_time(bgvalid,time_bg(nf),a14_time(10:14)
     1                    ,istatus)       

!      We might also check whether the file already resides on disk
       do iw = 1,n_written
         if(a14_time .eq. c_ftimes_written(iw))then
             write(6,*)' NOTE: skipping already processed file time '
     1                ,a14_time
             lga_status = 1
             goto900
         endif
       enddo ! iw

c Removal of this if/loop causes already existing lga files to be overwritten
c possibly with the same data.  However the error frequency on SBN may warrent
c this extra work.  
       if(nlga .gt. 0 .and. .false.)then
        do i=1,nlga
         print *,i,lga_names(i),':',fname_bg(nf),':',af_bg(nf),' test'       
         if (fname_bg(nf) .eq. lga_names(i)(bglen+1:bglen+13)) then
                  
           print *,i,lga_names(i),':',fname_bg(nf),':',af_bg(nf)
           call get_lgb_source(nx_laps,ny_laps             
     +                 ,fname_bg(nf),af_bg(nf),comment(1))

           if(cmodel(1:ic) .eq. comment(1)(1:ic)) then
              print *,'Note that LGB file exists. ',lga_names(i)
c             print *,'LGA file exists, not making one. ',lga_names(i)
c             lga_status=1
c             goto900
           else
              print *,'Overwritting LGA/LGB files ',lga_names(i)
     +             ,'from model ',comment(1)(1:i),' with new '
     +             ,'ones from model ',cmodel(1:ic)
           endif
         endif
        enddo          
       endif     

       if(bgpath(bglen:bglen).eq.'/')then
          fullname = bgpath(1:bglen)//fname_bg(nf)
       else
          fullname = bgpath(1:bglen)//'/'//fname_bg(nf)
       endif

       call s_len(fullname,i)
 
       print*
       print*,'Reading - ',fullname(1:i)
       print*


       allocate (htbg(nx_bg,ny_bg,nzbg_ht))
       allocate (tpbg(nx_bg,ny_bg,nzbg_tp))
       allocate (shbg(nx_bg,ny_bg,nzbg_sh))
       allocate (uwbg(nx_bg,ny_bg,nzbg_uv))
       allocate (vwbg(nx_bg,ny_bg,nzbg_uv))
       allocate (wwbg(nx_bg,ny_bg,nzbg_ww))
       allocate (prbght(nx_bg,ny_bg,nzbg_ht))
       allocate (prbguv(nx_bg,ny_bg,nzbg_uv))
       allocate (prbgsh(nx_bg,ny_bg,nzbg_sh))
       allocate (prbgww(nx_bg,ny_bg,nzbg_ww))

       allocate (htbg_sfc(nx_bg,ny_bg))
       allocate (prbg_sfc(nx_bg,ny_bg))
       allocate (shbg_sfc(nx_bg,ny_bg))
       allocate (uwbg_sfc(nx_bg,ny_bg))
       allocate (vwbg_sfc(nx_bg,ny_bg))
       allocate (tdbg_sfc(nx_bg,ny_bg))
       allocate (tpbg_sfc(nx_bg,ny_bg))
       allocate (t_at_sfc(nx_bg,ny_bg))
       allocate (mslpbg(nx_bg,ny_bg))
       allocate (pcpbg(nx_bg,ny_bg))
       allocate (crefbg(nx_bg,ny_bg))
       allocate (tpwbg(nx_bg,ny_bg))
       allocate (cwatbg(nx_bg,ny_bg))
       allocate (swibg(nx_bg,ny_bg))

       t_at_sfc=missingflag

       call read_bgdata(nx_bg,ny_bg
     +    ,nzbg_ht,nzbg_tp,nzbg_sh,nzbg_uv,nzbg_ww
     +    ,'lapsb',bgpath,fname_bg(nf),af_bg(nf)
     +    ,fullname,cmodel,bgmodel
     +    ,prbght,prbgsh,prbguv,prbgww
     +    ,htbg,tpbg,uwbg,vwbg,shbg,wwbg
     +    ,htbg_sfc,prbg_sfc,shbg_sfc,tdbg_sfc,tpbg_sfc
     +    ,t_at_sfc,uwbg_sfc,vwbg_sfc,mslpbg,pcpbg,crefbg
     +    ,tpwbg,cwatbg,swibg
     +    ,istatus_prep(nf))

       istatus=ishow_timer()
       write(6,*)' Returned from read_bgdata'

       if(bgmodel .eq. 13)then ! determine actual number of levels
           call get_plvls(plvl_grib, 100, nlvl_grib)
           write(6,*)' grib plvls info: ',nlvl_grib
     1                                   ,plvl_grib(1:nlvl_grib)
           nzbg_ht = 0
           do k = 1,nlvl_grib
               if(plvl_grib(k) .lt. 150000)then
                   nzbg_ht = nzbg_ht + 1
               endif
           enddo ! k

           write(6,*)' Reset nzbg_ht, etc. to ',nzbg_ht

           nzbg_tp=nzbg_ht
           nzbg_sh=nzbg_ht
           nzbg_uv=nzbg_ht
           nzbg_ww=nzbg_ht
       endif

       if(.false.)then
           print*,'After read'
           do k=1,nzbg_ww
              rmaxvv=maxval(wwbg(:,:,k))
              rminvv=minval(wwbg(:,:,k))
              print*,'k Max/Min vv ',k,rmaxvv,rminvv
           enddo
       endif

       if (istatus_prep(nf) .ne. 0) then

c           call s_len(fname_bg(nf),lf)
c           call s_len(bgpath,l)
c           if (bgmodel .gt. 1 .and. bgmodel .le. 3) then
c              fname13=fname_bg(nf)(1:lf)//af_bg(nf)
c           elseif (bgmodel .eq. 4) then
c              fname13=fname9_to_wfo_fname13(fname_bg(nf))
c           endif
          print *,'Background model data not returned from ',
     .'read_bgdata: ',bgpath(1:bglen)//fname_bg(nf)
          print *,'Process aborted for this file.'

c         convert to wfo if necessary

          print*
          print*,'Updating reject information: '

          if(bgmodel.eq.4 .or. (bgmodel.eq.10))then   !WNI -- added bgmodel 10
             print*,'WFO Update'
             fname13=wfo_fname13_to_fname9(fname_bg(nf))
             fname13=fname13(1:9)//af_bg(nf)
             do ii=1,accepted_files
                if(fname13.eq.bg_names(ii))then
                   reject_cnt=reject_cnt+1
                   reject_names(reject_cnt)=bg_names(ii)
                   print*,'reject_cnt/reject_names'
                   print*,'cnt/time: ',reject_cnt
     +,reject_names(reject_cnt)

                endif
             enddo
          else
             print*
             do ii=1,accepted_files
                if(fname_bg(nf).eq.bg_names(ii))then
                   reject_cnt=reject_cnt+1
                   reject_names(reject_cnt)=bg_names(ii)

                   print*,'reject_cnt/reject_names'
                   print*,'cnt/time: ',reject_cnt
     +,reject_names(reject_cnt)
                   print*,'fullname',fullname
                endif
             enddo
          endif

          lga_status= -nf

          deallocate(htbg, 
     +               tpbg, 
     +               shbg, 
     +               uwbg, 
     +               vwbg, 
     +               wwbg,
     +               prbght, 
     +               prbguv, 
     +               prbgsh, 
     +               prbgww, 
     +               htbg_sfc, 
!    +               prbg_sfc,
     +               shbg_sfc, 
     +               uwbg_sfc, 
     +               vwbg_sfc, 
     +               tdbg_sfc,
     +               tpbg_sfc,
     +               t_at_sfc,
     +               mslpbg,
     +               pcpbg,crefbg,tpwbg,cwatbg,swibg)

 
       else   !processing the file because it is not a reject

c
c ****** Vertically interpolate background data to LAPS isobaric levels.
c
         if(linterp)then   ! this switch determines if we are going to h/v-interp or not

           itstatus(1)=ishow_timer()

c ****** Horizontally interpolate background data to LAPS grid points. ********
c
!          itstatus(2)=ishow_timer()

!          Get bounding box of hinterp to restrict vinterp range
           write(6,*)' calling init_hinterp'
           call init_hinterp(nx_bg,ny_bg,nx_laps,ny_laps,gproj,
     .        lat,lon,grx,gry,bgmodel,cmodel,wrapped)

           print*,'LAPS (Input) Grid Corners'
           print*, 'SW: grx(1,1)/gry(1,1) ', grx(1,1),gry(1,1)
           print*, 'SE: grx(nx,1)/gry(nx,1) '
     +, grx(nx_laps,1),gry(nx_laps,1)
           print*, 'NW: grx(1,ny/gry(1,ny) ',
     + grx(1,ny_laps),gry(1,ny_laps)
           print*, 'NE: grx(nx,ny)/gry(nx,ny) '
     +, grx(nx_laps,ny_laps),gry(nx_laps,ny_laps)
           print*

c
c ***** Check if bkgd model domain satisfies analysis domain **********
c       ----------------------------------------------------
           lgrid_missing=.false.
           xmin = 999999.
           ymin = 999999.
           xmax =-999999.
           ymax =-999999.
           search_grid_missing: do i=1,nx_laps
            do j=1,ny_laps
             if(grx(i,j).eq.missingflag .or.
     .           gry(i,j).eq.missingflag)then
                 lgrid_missing=.true.
                 exit search_grid_missing
             else ! xygrid exists at LAPS grid point
                 xmin = min(grx(i,j),xmin)
                 xmax = max(grx(i,j),xmax)
                 ymin = min(gry(i,j),ymin)
                 ymax = max(gry(i,j),ymax)
             endif
            enddo
           enddo search_grid_missing

           print*,'LAPS (Input) Grid Bounding Box'
           print*,'Xrange ',xmin,xmax
           print*,'Yrange ',ymin,ymax

!          Add a perimeter to account for hinterp spline
           call get_grid_spacing_cen(grid_spacing_cen_m,istatus)
           bgres = 40000. ! worst case for now
!          i_perim = int(bgres / grid_spacing_cen_m) + 1
           i_perim = 5
           ixmin = max(nint(xmin)-i_perim,1)
           ixmax = min(nint(xmax)+i_perim,nx_bg)
           iymin = max(nint(ymin)-i_perim,1)
           iymax = min(nint(ymax)+i_perim,ny_bg)
           write(6,*)' bgres / i_perim ',bgres,i_perim

           print*
           print*,'LAPS (Input) Grid Bounding Box For Vinterp'
           print*,'Xrange ',ixmin,ixmax
           print*,'Yrange ',iymin,iymax
c
           if(ixmin .gt. ixmax .OR. iymin .gt. iymax)then
               write(6,*)' WARNING: No valid bounding box'
               lga_status = -nf
               goto 999 ! deallocate/return
           endif

           istatus=ishow_timer()

           allocate( htvi(nx_bg,ny_bg,nz_laps),   !Height (m)
     .               tpvi(nx_bg,ny_bg,nz_laps),   !Temperature (K)
     .               shvi(nx_bg,ny_bg,nz_laps),   !Specific humidity (kg/kg)
     .               uwvi(nx_bg,ny_bg,nz_laps),   !U-wind (m/s)
     .               vwvi(nx_bg,ny_bg,nz_laps),   !V-wind (m/s)
     .               wwvi(nx_bg,ny_bg,nz_laps))   !W-wind (pa/s)

           if(.true.)then

!            The 'htbg_sfc' model terrain can be used or calculated (if 
!            unavailable).

             write(6,*)' htbg_sfc range: ',minval(htbg_sfc)
     1                                    ,maxval(htbg_sfc)       

             if(minval(htbg_sfc) .eq. missingflag .OR.
     1          maxval(htbg_sfc) .eq. missingflag      )then
                 write(6,*)
     1   ' WARNING: htbg_sfc has missing values, check BG model terrain'
 
                 write(6,*)' prbg_sfc range: ',minval(prbg_sfc)
     1                                        ,maxval(prbg_sfc)       

                 if(.false.)then
                     write(6,*)
     1   ' Calculating model terrain from model pressure and std atmos'
                     do j = 1,ny_bg
                     do i = 1,nx_bg
                        htbg_sfc(i,j) = PsaToZ(prbg_sfc(i,j)/100.)
                     enddo ! i
                     enddo ! j                 
                 else
                     write(6,*)
     1   ' Calculating model terrain from model sfc pres and 3D height' 
                     call pres_to_ht_2d(prbg_sfc,prbght*100.,htbg,tpbg
     1                                 ,nx_bg,ny_bg,nzbg_ht,htbg_sfc
     1                                 ,istatus)
                 endif

                 write(6,*)' recalculated htbg_sfc range: '
     1                                        ,minval(htbg_sfc)
     1                                        ,maxval(htbg_sfc)       

             endif

           endif ! .true.

           if(vertical_grid .eq. 'PRESSURE')then 
       
             write(6,*)
             do k = 1,nzbg_ht
                 write(6,*)' htbg range at level ',k,minval(htbg(:,:,k))
     1                                              ,maxval(htbg(:,:,k))   
             enddo ! k

             write(6,*)
             do k = 1,nzbg_tp
                 write(6,*)' tpbg range at level ',k,minval(tpbg(:,:,k))
     1                                              ,maxval(tpbg(:,:,k))   
             enddo ! k

             call vinterp(nz_laps,nx_bg,ny_bg,nx_pr,ny_pr
     .         ,ixmin,ixmax,iymin,iymax
     .         ,nzbg_ht,nzbg_tp,nzbg_sh,nzbg_uv,nzbg_ww
     .         ,pr1d_mb,prbght,prbgsh,prbguv,prbgww
     .         ,htbg,tpbg,shbg,uwbg,vwbg,wwbg
     .         ,htvi,tpvi,shvi,uwvi,vwvi,wwvi)
       
             write(6,*)
             write(6,*)' compute tpvi range within LAPS box'
             do k = 1,nz_laps
                 write(6,*)' tpvi range at level ',k
     1                     ,minval(tpvi(ixmin:ixmax,iymin:iymax,k))
     1                     ,maxval(tpvi(ixmin:ixmax,iymin:iymax,k))   
             enddo ! k
       
             write(6,*)
             write(6,*)' compute htvi range within LAPS box'
             do k = 1,nz_laps
                 write(6,*)' htvi range at level ',k
     1                     ,minval(htvi(ixmin:ixmax,iymin:iymax,k))
     1                     ,maxval(htvi(ixmin:ixmax,iymin:iymax,k))   
             enddo ! k

           elseif(vertical_grid .eq. 'SIGMA_P')then
!            LAPS pressure should be on model grid (prvi grid). We want to 
!            model 'sigma_p' vertical levels on the model grid to construct 
!            the 3D 'prvi' array. The 'prbg_sfc' model terrain is used. 

             allocate(prvi(nx_bg,ny_bg,nz_laps))

             prtop = 100. ! Top of sigma grid in mb

             call get_sigmap_3d(prbg_sfc/100.,prtop,sigma1d,prvi
     1                         ,nx_bg,ny_bg,nz_laps,istatus)

             call vinterp(nz_laps,nx_bg,ny_bg,nx_pr,ny_pr
     .         ,ixmin,ixmax,iymin,iymax
     .         ,nzbg_ht,nzbg_tp,nzbg_sh,nzbg_uv,nzbg_ww
     .         ,prvi,prbght,prbgsh,prbguv,prbgww
     .         ,htbg,tpbg,shbg,uwbg,vwbg,wwbg
     .         ,htvi,tpvi,shvi,uwvi,vwvi,wwvi)

           elseif(vertical_grid .eq. 'SIGMA_HT')then
             allocate(prvi(nx_bg,ny_bg,nz_laps))

!            We want to model 'sigma_ht' vertical levels on the model 
!            horizontal grid to construct the 'htvi' array. 

             call get_ht_3d(nx_bg,ny_bg,nz_laps,htbg_sfc,htvi
     1                     ,istatus)
             if(istatus .ne. 1)then
                 write(6,*)' Error returned from get_ht_3d (model grid)'       
                 return
             endif
       
             do k = 1,nz_laps
                 write(6,*)' htvi range at level ',k,minval(htvi(:,:,k))
     1                                              ,maxval(htvi(:,:,k))   
             enddo ! k

             call get_ht_1d(nz_laps,ht_1d,istatus)
             if(istatus .ne. 1)then
                 write(6,*)' Error returned from get_ht_1d'
                 return
             endif

             write(6,*)
             do k = 1,nzbg_ht
                 write(6,*)' prbght range at level ',k
     1                                            ,minval(prbght(:,:,k))       
     1                                            ,maxval(prbght(:,:,k))   
             enddo ! i

             write(6,*)
             do k = 1,nzbg_ht
                 write(6,*)' htbg range at level ',k,minval(htbg(:,:,k))
     1                                              ,maxval(htbg(:,:,k))   
             enddo ! i

             write(6,*)
             do k = 1,nzbg_ht
                 write(6,*)' shbg range at level ',k,minval(shbg(:,:,k))
     1                                              ,maxval(shbg(:,:,k))   
             enddo ! i

             call vinterp_ht(nz_laps,nx_bg,ny_bg
     .         ,nzbg_ht,nzbg_tp,nzbg_sh,nzbg_uv,nzbg_ww
     .         ,htvi,prbght,prbgsh,prbguv,prbgww 
     .         ,htbg,tpbg,shbg,uwbg,vwbg,wwbg
     .         ,prvi,tpvi,shvi,uwvi,vwvi,wwvi) 

             write(6,*)
             do k = 1,nz_laps
                 write(6,*)' prvi range at level ',k,minval(prvi(:,:,k))       
     1                                              ,maxval(prvi(:,:,k))   
             enddo ! k

             write(6,*)
             do k = 1,nz_laps
                 write(6,*)' shvi range at level ',k,minval(shvi(:,:,k))       
     1                                              ,maxval(shvi(:,:,k))   
             enddo ! k

           else
             write(6,*)' Vertical grid not supported'
             stop

           endif

           itstatus(2)=ishow_timer()
           print*,' Vinterp elapsed time (sec): ',itstatus(2)
     1                                           -itstatus(1)

           deallocate (htbg, tpbg, shbg, uwbg, vwbg, wwbg
     +,prbght, prbguv, prbgsh, prbgww )

c
c ****** Run 2dx filter on vertically interpolated fields.
c ****** Only run filter on sh up to 300 mb, since the filter may
c           create small neg values when the field is small to begin with.
c ****** If bgmodel=4, there exists missing data.  Don't use these points
c           in the filter.
c
           allocate (msgpt(nx_bg,ny_bg))

           if (bgmodel .eq. 4 .or. bgmodel .eq. 9) then

            do j=1,ny_bg
            do i=1,nx_bg
               if (htvi(i,j,1) .eq. missingflag) then
                  msgpt(i,j)=0
               else
                  msgpt(i,j)=1
               endif
            enddo
            enddo
            do k=1,nz_laps
            do j=2,ny_bg-1
            do i=2,nx_bg-1
               if (max(msgpt(i-1,j-1),
     .                 msgpt(i-1,j  ),
     .                 msgpt(i-1,j+1),
     .                 msgpt(i  ,j-1),
     .                 msgpt(i  ,j  ),
     .                 msgpt(i  ,j+1),
     .                 msgpt(i+1,j-1),
     .                 msgpt(i+1,j  ),
     .                 msgpt(i+1,j+1)) .eq. 1) then
                  ct=0
                  htave=0.
                  tpave=0.
                  shave=0.
                  uwave=0.
                  vwave=0.
                  do jj=j-1,j+1
                  do ii=i-1,i+1
                     if (msgpt(ii,jj) .eq. 1) then
                        ct=ct+1
                        htave=htave+htvi(ii,jj,k)
                        tpave=tpave+tpvi(ii,jj,k)
                        shave=shave+shvi(ii,jj,k)
                        uwave=uwave+uwvi(ii,jj,k)
                        vwave=vwave+vwvi(ii,jj,k)
                     endif
                  enddo
                  enddo
                  cti=1./float(ct)
                  htvi(i,j,k)=htave*cti
                  tpvi(i,j,k)=tpave*cti
                  shvi(i,j,k)=shave*cti
                  uwvi(i,j,k)=uwave*cti
                  vwvi(i,j,k)=vwave*cti
               endif
            enddo
            enddo
            enddo
 
           endif ! bgmodel = 4 or 9
c
           do kk=nz_laps,1,-1
            if (pr1d_mb(kk) .ge. 300.) goto 20
           enddo
20         continue
           if (smooth_fields) THEN
             print*, "SMOOTHING BACKGROUND DATA!"
             call filter_2dx(htvi,nx_bg,ny_bg,nz_laps, 0.5)
             call filter_2dx(htvi,nx_bg,ny_bg,nz_laps,-0.5)
             call filter_2dx(tpvi,nx_bg,ny_bg,nz_laps, 0.5)
             call filter_2dx(tpvi,nx_bg,ny_bg,nz_laps,-0.5)
             call filter_2dx(shvi,nx_bg,ny_bg,kk     , 0.5)
             call filter_2dx(shvi,nx_bg,ny_bg,kk     ,-0.5)
             call filter_2dx(uwvi,nx_bg,ny_bg,nz_laps, 0.5)
             call filter_2dx(uwvi,nx_bg,ny_bg,nz_laps,-0.5)
             call filter_2dx(vwvi,nx_bg,ny_bg,nz_laps, 0.5)
             call filter_2dx(vwvi,nx_bg,ny_bg,nz_laps,-0.5)
             call filter_2dx(wwvi,nx_bg,ny_bg,nz_laps, 0.5)
             call filter_2dx(wwvi,nx_bg,ny_bg,nz_laps,-0.5)
           endif
c
           if (bgmodel .eq. 4) then
            do j=1,ny_bg
            do i=1,nx_bg
               if (msgpt(i,j) .eq. 0) then
                  do k=1,nz_laps
                     htvi(i,j,k)=missingflag
                     tpvi(i,j,k)=missingflag
                     shvi(i,j,k)=missingflag
                     uwvi(i,j,k)=missingflag
                     vwvi(i,j,k)=missingflag
                  enddo
               endif
            enddo
            enddo
           endif
c
c ***** Check if t_at_sfc is defined **********
c       ----------------------------
           lhif_tsfc=.true.
           search_tsfc_missing: do i=1,nx_bg
            do j=1,ny_bg
             if(t_at_sfc(i,j).eq.missingflag)then
                lhif_tsfc=.false.
                exit search_tsfc_missing
             endif
            enddo
           enddo search_tsfc_missing

c
          if(lgrid_missing)then

           print*,'Error: bkgd domain size insufficient for: '
     1           ,trim(cmodel)
           goto 999 ! return

          else
                 
!          Obtain height and pressure fields as needed on LAPS grid
           if(vertical_grid .ne. 'SIGMA_HT')then ! PRESSURE or SIGMA_P
              itstatus(2)=ishow_timer()
              print*,'use hinterp_field_3d for HT ',cmodel(1:ic)
              call hinterp_field_3d(nx_bg,ny_bg,nx_laps,ny_laps,nz_laps,       
     .                           grx,gry,htvi,ht,wrapped)

              write(6,*)
              do k = 1,nz_laps
                 write(6,*)' ht range at level ',k,minval(ht(:,:,k))       
     1                                            ,maxval(ht(:,:,k))   
              enddo ! k

              if(minval(ht) .eq. missingflag .OR.
     1           maxval(ht) .eq. missingflag      )then
                 write(6,*)
     1               ' ERROR: missing value in interpolated ht field'
              endif

              if(vertical_grid .eq. 'SIGMA_P')then ! generate 3D P     
                 call hinterp_field_3d(nx_bg,ny_bg,nx_laps,ny_laps,1,
     .                              grx,gry,prbg_sfc,pr_sfc,wrapped)

!                Convert 1D sigma levels to 3D pressures (with surface pressure)
                 call get_sigmap_3d(pr_sfc,prtop,sigma1d,prgd_pa
     1                             ,nx_laps,ny_laps,nz_laps,istatus)       
              endif

           else ! SIGMA_HT
              call get_ht_3d(nx_laps,ny_laps,nz_laps,topo,ht
     1                      ,istatus)
              if(istatus .ne. 1)then
                 write(6,*)' Error returned from get_ht_3d (LAPS grid)'
                 return
              endif

!             Input pressure (prvi) is being converted from mb to Pa
              call hinterp_field_3d(nx_bg,ny_bg,nx_laps,ny_laps,nz_laps, 
     .                           grx,gry,prvi*100.,prgd_pa,wrapped)

           endif

           if(.not. l_bilinear) then                                  
              call hinterp_field_3d(nx_bg,ny_bg,nx_laps,ny_laps,nz_laps,
     .                           grx,gry,uwvi,uw,wrapped)
              call hinterp_field_3d(nx_bg,ny_bg,nx_laps,ny_laps,nz_laps,
     .                           grx,gry,vwvi,vw,wrapped)
              call hinterp_field_3d(nx_bg,ny_bg,nx_laps,ny_laps,nz_laps,
     .                           grx,gry,tpvi,tp,wrapped)
              call hinterp_field_3d(nx_bg,ny_bg,nx_laps,ny_laps,nz_laps,
     .                           grx,gry,shvi,sh,wrapped)

              write(6,*)
              do k = 1,nz_laps
                 write(6,*)' sh range at level ',k,minval(sh(:,:,k))       
     1                                            ,maxval(sh(:,:,k))   
              enddo ! k

           else
              istatus=ishow_timer()
              print*,'use bilinear_laps_3d for U ',cmodel(1:ic)
              call bilinear_laps_3df(grx,gry,nx_bg,ny_bg
     .                             ,nx_laps,ny_laps,nz_laps,uwvi,uw)

              istatus=ishow_timer()
              print*,'use bilinear_laps_3d for V ',cmodel(1:ic)
              call bilinear_laps_3df(grx,gry,nx_bg,ny_bg
     .                             ,nx_laps,ny_laps,nz_laps,vwvi,vw)

              istatus=ishow_timer()
              print*,'use bilinear_laps_3d for T ',cmodel(1:ic)
              call bilinear_laps_3df(grx,gry,nx_bg,ny_bg
     .                             ,nx_laps,ny_laps,nz_laps,tpvi,tp)

              istatus=ishow_timer()        
              print*,'use bilinear_laps_3d for Q ',cmodel(1:ic)
              call bilinear_laps_3df(grx,gry,nx_bg,ny_bg  
     .                             ,nx_laps,ny_laps,nz_laps,shvi,sh)
           endif

           if(.not. lgb_only)then ! skip ww for sfc only option
                call hinterp_field_3d(nx_bg,ny_bg,
     .             nx_laps,ny_laps,nz_laps,grx,gry,wwvi,ww,wrapped)
           endif
          
           itstatus(3)=ishow_timer()
           print*,'Hinterp (3D) elapsed time (sec): ',itstatus(3)
     1                                               -itstatus(2)
c
c ****** Check for missing value flag in any of the fields.
c ****** Check for NaN's in any of the fields.
c
           deallocate(htvi,   !Height (m)
     .                tpvi,   !Temperature (K)
     .                shvi,   !Specific humidity (kg/kg)
     .                uwvi,   !U-wind (m/s)
     .                vwvi,   !V-wind (m/s)
     .                wwvi,   !W-wind (pa/s)
     .                msgpt)
        
           if(vertical_grid .eq. 'SIGMA_HT' .or. 
     1        vertical_grid .eq. 'SIGMA_P'       )then
               deallocate(prvi)   !3-D Interpolated Pressure (mb)
           endif

           do k=1,nz_laps
            do j=1,ny_laps
               do i=1,nx_laps
                  if((abs(ht(i,j,k)) .gt. 100000.) .or.
     +                   (ht(i,j,k)  .lt.-3000)  .or.
     +                   (tp(i,j,k)  .gt. 500.) .or.
     +                   (tp(i,j,k)  .lt. 150.) .or.
     +               (abs(sh(i,j,k)) .gt. 1.) .or.
     +               (abs(uw(i,j,k)) .gt. 150.) .or.
     +               (abs(vw(i,j,k)) .gt. 150.)       )then

c    +             .or.(abs(ww(i,j,k)) .gt. 10.)      )then
c
c ww may be missing from some models! Don't stop just because of that.
c                  if (max(ht(i,j,k),tp(i,j,k),sh(i,j,k),
c    .                     uw(i,j,k),vw(i,j,k)) .ge. missingflag)then       

                     print*,
     + 'ERROR: interpolated values exceed allowed ranges at gridpoint: '
     +                      ,i,j,k           
                     print*,'ht/tp/sh/uw/vw/ww: ',ht(i,j,k),tp(i,j,k)        
     +                      ,sh(i,j,k), uw(i,j,k),vw(i,j,k),ww(i,j,k)

                     if(bgmodel .eq. 13)then
                        print*,'Check input model data with GRIB viewer'
                     else
                        print*,'Check input model data'
                     endif

                     reject_cnt=reject_cnt+1
!                    reject_names(reject_cnt)=bg_names(ii)

                     print*,'reject_cnt/reject_names'
                     print*,'cnt/time: ',reject_cnt
     +                     ,' UNKNOWN reject_name'
!    +                     ,reject_names(reject_cnt)

                     lga_status = -nf

                     deallocate (htbg_sfc
     +                          ,prbg_sfc
     +                          ,shbg_sfc
     +                          ,uwbg_sfc
     +                          ,vwbg_sfc
     +                          ,tdbg_sfc
     +                          ,tpbg_sfc
     +                          ,t_at_sfc
     +                          ,mslpbg
     +                          ,pcpbg,crefbg,tpwbg,cwatbg,swibg)

                     goto 999 ! deallocate/return

                  endif
               enddo
            enddo
           enddo
c
           call check_nan3 (ht,nx_laps,ny_laps,nz_laps,nan_flag)
           if(nan_flag .ne. 1) then
            print *,' ERROR: NaN found in array ht'
            lga_status = -nf
            goto 999 ! deallocate/return
           endif
c
           call check_nan3 (tp,nx_laps,ny_laps,nz_laps,nan_flag)
           if(nan_flag .ne. 1) then
            print *,' ERROR: NaN found in array tp'
            lga_status = -nf
            goto 999 ! deallocate/return
           endif
c
           call check_nan3 (sh,nx_laps,ny_laps,nz_laps,nan_flag)
           if(nan_flag .ne. 1) then
            print *,' ERROR: NaN found in array sh'
            lga_status = -nf
            goto 999 ! deallocate/return
           endif
c
           call check_nan3 (uw,nx_laps,ny_laps,nz_laps,nan_flag)
           if(nan_flag .ne. 1) then
            print *,' ERROR: NaN found in array uw'
            lga_status = -nf
            goto 999 ! deallocate/return
           endif
c
           call check_nan3 (vw,nx_laps,ny_laps,nz_laps,nan_flag)
           if(nan_flag .ne. 1) then
            print *,' ERROR: NaN found in array vw'
            lga_status = -nf
            goto 999 ! deallocate/return
           endif
c
           call check_nan3 (ww,nx_laps,ny_laps,nz_laps,nan_flag)
           if(nan_flag .ne. 1) then
            print *,' ERROR: NaN found in array vw'
            lga_status = -nf
            goto 999 ! deallocate/return
           endif
c
c ****** Horizontally interpolate background surface data to LAPS grid points.
c
           itstatus(3)=ishow_timer()

!          if(bgmodel.ne.1.and.bgmodel.ne.9 .AND.
!    1                         bgmodel.ne.3      )then
           if(.true.)then

            write(6,*)' Calling hinterp_field for surface variables'

            if(minval(htbg_sfc) .eq. missingflag .OR.
     1         maxval(htbg_sfc) .eq. missingflag      )then
              write(6,*)' ERROR: htbg_sfc has missing data'
              return
            else
              call hinterp_field_2d(nx_bg,ny_bg,nx_laps,ny_laps,1,
     .        grx,gry,htbg_sfc,ht_sfc,wrapped)
            endif

            if(minval(tdbg_sfc) .eq. missingflag .OR.
     1         maxval(tdbg_sfc) .eq. missingflag      )then
              write(6,*)' NOTE: tdbg_sfc has missing data'
              td_sfc = missingflag
            else
              call hinterp_field_2d(nx_bg,ny_bg,nx_laps,ny_laps,1,
     .        grx,gry,tdbg_sfc,td_sfc,wrapped)
            endif

            call hinterp_field_2d(nx_bg,ny_bg,nx_laps,ny_laps,1,
     .        grx,gry,tpbg_sfc,tp_sfc,wrapped)

            if(minval(shbg_sfc) .eq. missingflag .OR.
     1         maxval(shbg_sfc) .eq. missingflag      )then
              write(6,*)' NOTE: shbg_sfc has missing data'
              sh_sfc = missingflag
            else
              call hinterp_field_2d(nx_bg,ny_bg,nx_laps,ny_laps,1,
     .        grx,gry,shbg_sfc,sh_sfc,wrapped)
            endif

            call hinterp_field_2d(nx_bg,ny_bg,nx_laps,ny_laps,1,
     .        grx,gry,uwbg_sfc,uw_sfc,wrapped)
            call hinterp_field_2d(nx_bg,ny_bg,nx_laps,ny_laps,1,
     .        grx,gry,vwbg_sfc,vw_sfc,wrapped)

            if(minval(prbg_sfc) .eq. missingflag .OR.
     1         maxval(prbg_sfc) .eq. missingflag      )then
              write(6,*)' WARNING: prbg_sfc has missing data'
              pr_sfc = missingflag
            else
              call hinterp_field_2d(nx_bg,ny_bg,nx_laps,ny_laps,1,
     .        grx,gry,prbg_sfc,pr_sfc,wrapped)
            endif

            call hinterp_field_2d(nx_bg,ny_bg,nx_laps,ny_laps,1,
     .        grx,gry,mslpbg,mslp,wrapped)
            call hinterp_field_2d(nx_bg,ny_bg,nx_laps,ny_laps,1,
     .        grx,gry,pcpbg,pcp_sfc,wrapped)
            call hinterp_field_2d(nx_bg,ny_bg,nx_laps,ny_laps,1,
     .        grx,gry,cwatbg,cw_sfc,wrapped)

            write(6,*)' ht_sfc range = ',minval(ht_sfc),maxval(ht_sfc)
            write(6,*)' pr_sfc range = ',minval(pr_sfc),maxval(pr_sfc)
            write(6,*)' td_sfc range = ',minval(td_sfc),maxval(td_sfc)
            write(6,*)' sh_sfc range = ',minval(sh_sfc),maxval(sh_sfc)
            write(6,*)' mslp range = ',minval(mslp),maxval(mslp)
            write(6,*)' cw_sfc range = ',minval(cw_sfc),maxval(cw_sfc)

            if(maxval(mslp) .le. 70000.) then
              print *,' ERROR: MSLP out of allowed range'
              lga_status = -nf
              goto 999 ! deallocate/return
            endif

            if(maxval(pr_sfc) .le. 30000.) then
              print *,' ERROR: pr_sfc out of allowed range'
              lga_status = -nf
              goto 999 ! deallocate/return
            endif
c
c Because not all model backgrounds have t_at_sfc (ground and/or sst)
c then no need to hinterp unless it exists.
c
            if(lhif_tsfc)then
               print*,'Horizontally Interpolate T at Sfc (tgd)'        
               call hinterp_field_2d(nx_bg,ny_bg,nx_laps,ny_laps,1,
     .           grx,gry,t_at_sfc,t_sfc,wrapped)
            else
               print*,'DO NOT Horizontally Interpolate T at Sfc (tgd)'       
               t_sfc=missingflag
            endif

            itstatus(3)=ishow_timer()
c
c check for T > Td before sfc p computation. Due to large scale
c interpolation we can have slightly larger (fractional) Td than T.
c
            call tdcheck(nx_laps,ny_laps,td_sfc,tp_sfc,
     &icnt,i_mx,j_mx,i_mn,j_mn,diff_mx,diff_mn)

            print *,' Dewpoint check (before call sfcbkgd):'
            print *,'     Dewpt greater than temp at ',icnt,' points.'

            if(icnt .gt. 0) then
               print*,' Max diff of ',diff_mx,' at ',i_mx,',',j_mx
               print*,' Min diff of ',diff_mn,' at ',i_mn,',',j_mn
            endif

            write(6,*)' td_sfc range = ',minval(td_sfc),maxval(td_sfc)

            itstatus(3)=ishow_timer()

           else ! bgmodel = 1,3,9

            call hinterp_field_2d(nx_bg,ny_bg,nx_laps,ny_laps,1,
     .                  grx,gry,prbg_sfc,pr_sfc,wrapped)
            call hinterp_field_2d(nx_bg,ny_bg,nx_laps,ny_laps,1,
     .                  grx,gry,mslpbg,mslp,wrapped)

           endif ! bgmodel .NE. 1,3,9

           deallocate (htbg_sfc)
           deallocate (prbg_sfc)
           deallocate (shbg_sfc)
           deallocate (uwbg_sfc)
           deallocate (vwbg_sfc)
           deallocate (tdbg_sfc)
           deallocate (tpbg_sfc)
           deallocate (t_at_sfc) 
           deallocate (mslpbg)
           deallocate (pcpbg,crefbg,tpwbg,cwatbg,swibg)
c
c... Do the temp, moisture (sh_sfc returns with Td), and pressures
c... Only ETA48_CONUS and namelist switch "luse_sfc_bkgd" enable the use
c... of subroutine sfcbkgd_sfc. This routine uses the 2m Td and sfc_press
c... 2D arrays directly from the background model.

           if(luse_sfc_bkgd)then ! tested only for ETA48_CONUS
             if(minval(tp_sfc) .eq. missingflag .AND.
     1          minval(tp_sfc) .eq. missingflag       )then
               write(6,*)' ERROR - surface temperature is missing'
               return
             endif
             if(vertical_grid .eq. 'PRESSURE')then 
               call sfcbkgd_sfc(bgmodel,tp,sh,ht,ht_sfc,td_sfc,tp_sfc
     .           ,sh_sfc,topo,pr1d_pa,nx_laps, ny_laps, nz_laps, pr_sfc
     .           ,nx_pr,ny_pr)      
               td_sfc_hi = td_sfc
             elseif(vertical_grid .eq. 'SIGMA_P' .or.
     .              vertical_grid .eq. 'SIGMA_HT'      )then
               call sfcbkgd_sfc(bgmodel,tp,sh,ht,ht_sfc,td_sfc,tp_sfc
     .           ,sh_sfc,topo,prgd_pa,nx_laps, ny_laps, nz_laps, pr_sfc
     .           ,nx_pr,ny_pr)      
               td_sfc_hi = td_sfc
             endif

           else
!             write(6,*)' td_sfc range = ',minval(td_sfc),maxval(td_sfc)
              if(vertical_grid .eq. 'SIGMA_HT')then
                write(6,*)
     1           ' WARNING: not using surface fields for SIGMA_HT grid'       
              endif
              call sfcbkgd(bgmodel,tp,sh,ht,tp_sfc,sh_sfc,td_sfc
     .           ,td_sfc_hi, topo
     .           ,pr1d_pa, nx_laps, ny_laps, nz_laps, pr_sfc
     .           ,nx_pr,ny_pr,istatus)      

           endif

           write(6,*)' td_sfc range = ',minval(td_sfc),maxval(td_sfc)

           call tdcheck(nx_laps,ny_laps,td_sfc,tp_sfc,
     &icnt,i_mx,j_mx,i_mn,j_mn,diff_mx,diff_mn)
           print *,' Dewpoint check (after call sfcbkgd):'
           print *,'     Dewpt greater than temp at ',icnt,' points.'

           if(icnt .gt. 0) then
            print*,'Max diff of ',diff_mx,' at ',i_mx,',',j_mx
            print*,'Min diff of ',diff_mn,' at ',i_mn,',',j_mn
c
c fix sfc Td to not be greater than T at points determined above
            where(td_sfc .gt. tp_sfc)td_sfc=tp_sfc
           endif
c
c..... Do the winds
c
           itstatus(3)=ishow_timer()
           if(ltest_vertical_grid('PRESSURE'))then
               write(6,*)' Interpolate 3D winds to the hi-res surface'
               call interp_to_sfc(topo,uw,ht,nx_laps,ny_laps,
     &                            nz_laps,missingflag,uw_sfc)
               call interp_to_sfc(topo,vw,ht,nx_laps,ny_laps,
     &                            nz_laps,missingflag,vw_sfc)
           elseif(ltest_vertical_grid('SIGMA_P'))then
               write(6,*)' Use lowest level 3-D sigma winds for the sfc'
               if(sigma1d(1) .eq. 1.0)then ! lowest sigma is at the sfc
                   uw_sfc(:,:) = uw(:,:,1)
                   vw_sfc(:,:) = vw(:,:,1)
               else
                   write(6,*)' ERROR: Unable to interpolate wind to sfc'
                   stop
               endif
           else ! SIGMA_HT
               write(6,*)' Use lowest level 3-D sigma winds for the sfc'
               if(.true.)then ! assume lowest sigma is at the sfc
                   uw_sfc(:,:) = uw(:,:,1)
                   vw_sfc(:,:) = vw(:,:,1)
               else
                   write(6,*)' ERROR: Unable to interpolate wind to sfc'
                   stop
               endif
           endif

           itstatus(3)=ishow_timer()

c
c..... Compute reduced pressure using reduced pressure level from
c      surface namelist file
c
c          Read surface parameters into module memory structure
           call get_directory('static',static_dir,len_dir)
           filename = static_dir(1:len_dir)//'/surface_analysis.nl'
c          call read_namelist_laps('sfc_anal',filename)
           call get_laps_redp(redp_lvl,istatus)

           allocate (rp_lvl(nx_laps,ny_laps)
     1              ,rp_tp(nx_laps,ny_laps),rp_sh(nx_laps,ny_laps)
     1              ,rp_td(nx_laps,ny_laps))

           do j=1,ny_laps
           do i=1,nx_laps
              rp_lvl(i,j)=redp_lvl
              rp_tp(i,j)=tp_sfc(i,j)
              rp_sh(i,j)=sh_sfc(i,j)                                 
              rp_td(i,j)=td_sfc(i,j)                                 
!             rp_sh(i,j)=td_sfc(i,j) ! borrow rp_sh for passing td_sfc
           enddo
           enddo
c
c always use sfcbkgd (as opposed to sfcbkgd_sfc) to compute reduced pressure
c because this version uses the 3D analysis info for computations.
c
           itstatus(3)=ishow_timer()

!          if(vertical_grid .ne. 'SIGMA_HT')then
           write(6,*)' call sfcbkgd for reduced pressure'
           if(vertical_grid .eq. 'PRESSURE')then 
              call sfcbkgd(bgmodel,tp,sh,ht,rp_tp,rp_sh
     1                    ,rp_td,dum2_2d,rp_lvl,pr1d_pa
     1                    ,nx_laps, ny_laps, nz_laps, rp_sfc
     1                    ,nx_pr,ny_pr,istatus)      
           else
               call sfcbkgd(bgmodel,tp,sh,ht,rp_tp,rp_sh
     1                    ,rp_td,dum2_2d,rp_lvl,prgd_pa
     1                    ,nx_laps, ny_laps, nz_laps, rp_sfc
     1                    ,nx_pr,ny_pr,istatus)      
           endif

           deallocate (rp_lvl,rp_tp,rp_sh,rp_td)

           itstatus(3)=ishow_timer()
c
c ****** Eliminate any supersaturations or negative sh generated 
c           through interpolation (set min sh to 1.e-6).
c
           write(6,*)' Eliminate supersaturations'
           icnt=0
           do k=1,nz_laps
           do j=1,ny_laps
           do i=1,nx_laps
            if(tp(i,j,k).gt.100.0)then
               if(vertical_grid .ne. 'SIGMA_HT')then
                 shsat=ssh2(pr1d_mb(k),tp(i,j,k)-273.15,
     .             tp(i,j,k)-273.15,-132.0)*0.001
               else
                 shsat=ssh2(prgd_pa(i,j,k)/100.,tp(i,j,k)-273.15,
     .             tp(i,j,k)-273.15,-132.0)*0.001
               endif
               sh(i,j,k)=max(1.0e-6,min(sh(i,j,k),shsat))
            else
               icnt=icnt+1
            endif
           enddo
           enddo
           enddo
           print*
           if(icnt.gt.0)then
            print*,'Warning: found ',icnt,' 3D points when'
            print*,'         checking for supersaturations'
           endif
c
           itstatus(4)=ishow_timer()
           print*,'Hinterp (2D) elapsed time (sec): ',itstatus(4)
     1                                               -itstatus(3)
           print*
c
c the wind components are still on the native grid projection;
c rotate them to the LAPS (output) domain as necessary.

!          itstatus_rot=ishow_timer()

           call rotate_background_uv(nx_laps,ny_laps,nz_laps,lon
     &,bgmodel,cmodel,fullname,gproj,lon0,lat0,lat1,uw,vw,uw_sfc,vw_sfc
     &,lgb_only,istatus)
           if(istatus.ne.1)then
              print*,'Error in rotate_background_uv '
              return
           endif

           istatus=ishow_timer()
!          print*,'After rotation: elapsed time (sec): ',itstatus_rot
!          print*

          endif !if grx/gry from init_hinterp are defined (lgrid_missing)

         else !linterp is false

c this is a grid compatible fua file

          ht=htbg
          tp=tpbg
          sh=shbg
          uw=uwbg
          vw=vwbg
          ww=wwbg
          pr_sfc=prbg_sfc
          mslp=mslpbg
c         td_sfc=tdbg_sfc
          td_sfc=shbg_sfc
          tp_sfc=tpbg_sfc
          ht_sfc=htbg_sfc
          sh_sfc=shbg_sfc
          uw_sfc=uwbg_sfc
          vw_sfc=vwbg_sfc
          pcp_sfc=pcpbg
          cw_sfc=cwatbg
c
c LAPS_FUA doesn't require interp but we still want to recompute
c pr_sfc, tp_sfc and sh_sfc using high res terrain
c

          if(cmodel.eq.'LAPS_FUA')then ! SIGMA_P option needed  
              call sfcbkgd_sfc(bgmodel,tp,sh,ht,ht_sfc
     &,td_sfc,tp_sfc,sh_sfc,topo,pr1d_mb,nx_laps,ny_laps,nz_laps,pr_sfc
     &                        ,nx_pr,ny_pr)      
              call tdcheck(nx_laps,ny_laps,td_sfc,tp_sfc,
     &icnt,i_mx,j_mx,i_mn,j_mn,diff_mx,diff_mn)

              print *,' Td check (after call sfcbkgd - LAPS_FUA):'
              print *,' Td greater than T at ',icnt,' points.'
              if(icnt .gt. 0) then
                 print*,'Max diff = ',diff_mx,' at ',i_mx,',',j_mx
                 print*,'Min diff = ',diff_mn,' at ',i_mn,',',j_mn
c
c fix sfc Td to not be greater than T at points determined above
c
                 where(td_sfc .gt. tp_sfc)td_sfc=tp_sfc
              endif
c
c..... Do the winds
c
              write(6,*)' Interpolate 3D winds to the hi-res surface'
              call interp_to_sfc(topo,uw,ht,nx_laps,ny_laps,
     &                         nz_laps,missingflag,uw_sfc)
              call interp_to_sfc(topo,vw,ht,nx_laps,ny_laps,
     &                         nz_laps,missingflag,vw_sfc)

          endif !(cmodel .eq. 'LAPS_FUA')

          deallocate (htbg, tpbg, shbg, uwbg, vwbg, wwbg
     +        ,prbght, prbguv, prbgsh, prbgww )
          deallocate (htbg_sfc,prbg_sfc,shbg_sfc,uwbg_sfc
     +        ,vwbg_sfc,tdbg_sfc, tpbg_sfc, t_at_sfc, mslpbg, pcpbg
     +        ,crefbg, tpwbg, cwatbg, swibg)

         endif !(linterp)

         itstatus(5)=ishow_timer()
c
c Write LGA
c ---------
         write(6,*)' Writing lga and/or lgb: ',trim(cmodel),' '
     1                 ,time_bg(nf),valid_bg(nf),bgvalid
         if(.not.lgb_only)then

          if(vertical_grid .eq. 'PRESSURE')then
            call write_lga(nx_laps,ny_laps,nz_laps,time_bg(nf),
     .bgvalid,cmodel,missingflag,pr1d_mb,ht,tp,sh,uw,vw,ww,istatus)
            if(istatus.ne.1)then
             print*,'Error writing lga - returning to main'
             return
            endif
          elseif(vertical_grid .eq. 'SIGMA_P')then 
            call write_lga(nx_laps,ny_laps,nz_laps,time_bg(nf),
     .bgvalid,cmodel,missingflag,sigma1d,ht,tp,sh,uw,vw,ww,istatus)
            if(istatus.ne.1)then
             print*,'Error writing lga - returning to main'
             return
            endif
          elseif(vertical_grid .eq. 'SIGMA_HT')then 
            call write_lgap(nx_laps,ny_laps,nz_laps,time_bg(nf),
     .bgvalid,cmodel,missingflag,ht_1d,prgd_pa,tp,sh,uw,vw,ww,istatus)
            if(istatus.ne.1)then
             print*,'Error writing lga - returning to main'
             return
            endif
          endif

         endif
c         
c Write LGB
c ---------

         if(bgmodel.ne.7)then
          do j=1,ny_laps
          do i=1,nx_laps
            if(pr_sfc(i,j) .lt. missingflag) then
               td_sfc_c = td_sfc_hi(i,j)-273.15
               if(td_sfc_c .lt. -199.)then
                   write(6,*)' ERROR, td_sfc_c < -199.'
                   istatus = 0
                   return
               endif

               tp_sfc_c = tp_sfc(i,j)-273.15
               if(tp_sfc_c .lt. td_sfc_c)then
                   write(6,*)' WARNING: TSfc < TdSfc ',tp_sfc_c,td_sfc_c
               endif
               qsfc(i,j)=ssh2(pr_sfc(i,j)*0.01,
     +                   tp_sfc_c,
     +                   td_sfc_c,-132.)*0.001
c              sfcgrid(i,j,kk+4)=qsfc(i,j)

            else
               qsfc(i,j) = missingflag
c              sfcgrid(i,j,kk+4)=missingflag

            endif

          enddo
          enddo
          call write_lgb(nx_laps,ny_laps,time_bg(nf),bgvalid
     .,cmodel,missingflag,uw_sfc,vw_sfc,tp_sfc,t_sfc,qsfc
     .,pr_sfc,mslp,td_sfc_hi,rp_sfc,pcp_sfc,cw_sfc,istatus)
          if(istatus.ne.1)then
            print*,'Error writing lgb - returning to main'
            return
          endif

          itstatus(6)=ishow_timer()
          print*,'Elapsed time - write grids (sec): ',itstatus(6)
     1                                               -itstatus(5)
          print*

         endif

         n_written = n_written + 1
         c_ftimes_written(n_written) = a14_time

         lga_status = 1
c
       endif !istatus_prep(nf) -> if the background file was properly read
                                 !and no bad data was found

 900   continue

      enddo ! Main loop through two model backgrounds (nf)

      if (allocated(prgd_pa))deallocate (prgd_pa)

c time interpolate between existing lga (bg) files.
c-------------------------------------------------------------------------------
c
      if(lga_status.eq.1 .and. mode .gt. 1)then
       itstatus(7)=ishow_timer()
       do i=1,nbg
          call s_len(bg_names(i),j)
          print*,'bg_name: ',i,bg_names(i)(1:j)
       enddo
c
c *** Determine if new file needs to (can) be created and perform
c        linear time interpolation.
c
       if(accepted_files.gt.1)then

         i=accepted_files
c        if(accepted_files.gt.2)then
c           i=2
c        endif

         print*,i,bg_times(i),bg_times(i-1),
     +     bg_valid(i),bg_valid(i-1),laps_cycle_time,
     +     i4time_bg_valid(i),i4time_bg_valid(i-1)

!        i4time_bg_valid times expected to be in reverse order 
         if(i4time_bg_valid(i-1) .gt.i4time_now   .and.
     +      i4time_bg_valid(i)   .lt.i4time_now  )then
            ext = 'lga'
            call get_directory(ext,outdir,len_dir) 
            print*,outdir,ext,nz_laps

c interp 3D fields
            if(.not.lgb_only)then

               call time_interp(outdir,ext,
     +           nx_laps,ny_laps,nz_laps,6,pr1d_mb,ht_1d,
     +           i4time_bg_valid(i),i4time_bg_valid(i-1),
     +           i4time_now,bg_times(i-1),bg_valid(i-1),
     +           bg_times(i  ),bg_valid(i  ))

            endif

            if(bgmodel.ne.1.or.bgmodel.ne.7)then
               ext = 'lgb'
               call get_directory(ext,outdir,len_dir) 
               print*,outdir,ext

c interp 2D fields
               call time_interp(outdir,ext,
     +           nx_laps,ny_laps,1,10,pr1d_mb(1),ht_1d(1),
     +           i4time_bg_valid(i),i4time_bg_valid(i-1),
     +           i4time_now,bg_times(i-1),bg_valid(i-1),
     +           bg_times(i  ),bg_valid(i  ))
            endif
         elseif(cmodel.eq.'LAPS')then
c            if(i4time_bg_valid(i-1) .ge.i4time_now   .and.
c    +          i4time_bg_valid(i)   .le.i4time_now  )then
                print*,'Determine how to time interpolate'
                print*,'and advance anal for t+1 cycle'
c            endif
         else
            print*,'Time Interpolation Not Necessary!'
         endif

         itstatus(8)=ishow_timer()
         print*,'Elapsed time - interp grids (sec): ',itstatus(8)
     1                                               -itstatus(7)
         print*

       elseif(accepted_files.eq.1)then

c code isn't ready to go yet.
c        if(i4time_bg_valid(i).gt.i4time_now)then
c           do k=1,nlga
c              lgadiff=i4time_bg_valid(i)-lga_times(k)
c           enddo
c        endif 

         print*,'Determine if existing lga plus the new one ',
     .          'can be used for time interpolation'
         print*
         print*,'Software does not exist yet. You lose.'

       else
         print*,'No time interp when accepted_files = 0'
         print*
       endif

      elseif(lga_status.eq.0 .and. nbg.eq.0)then
       print*,'bkgd files already exists; no new ones to proc'
c      lga_status = 1
      elseif(mode .le. 1)then
       print*,'spatial interp only was selected, skip time interp'       
      else
       print*,'No time interp when bad data found'
       print*
      endif

      print*,'Total time elapsed lga: '
!     print*,'(sec) :',itstatus(1)+itstatus(2)+itstatus(3)
!    &+itstatus(4)+itstatus(5)
      itstatus(10) = ishow_timer()

 999  if (allocated(prgd_pa)) deallocate(prgd_pa)
      if (allocated(prbg_sfc))deallocate(prbg_sfc)
      if (allocated(prvi))deallocate(prvi)

      return
      end
c
c=========================================================
c
      subroutine get_lga_source(nx,ny,nz,fname,af,source)
      character*9 fname
      character*4 af
      character*10 dumb1
      character*125 comment
      character*12 source
      integer ihour, bgtime, nx,ny,nz
      real dumb2(nx,ny,nz)

      call i4time_fname_lp(fname,bgtime,istatus)
      read(af,'(i4)') ihour
      call get_lapsdata_3d(bgtime,bgtime+ihour*3600,nx,ny,
     +        nz,'lga ','HT ',dumb1,comment,dumb2,istatus)
      if(istatus.lt.1) then
         stop 'error returned from get_lapsdata_3d'
      endif


      source = comment(1:12)
      return
      end
c
c=========================================================
c
      subroutine get_lgb_source(nx,ny,fname,af,source)
      character*9 fname
      character*4 af
      character*10 units
      character*125 comment
      character*12 source
      integer ihour, bgtime, nx,ny
      real dumb2(nx,ny)

      call i4time_fname_lp(fname,bgtime,istatus)
      read(af,'(2i2)') ihour,imin
      bgvalidtime = bgtime + ihour*3600 + imin*60
      call get_lapsdata_2d(bgtime,bgvalidtime,'lgb ','TSF',units
     .                    ,comment,nx,ny,dumb2,istatus)
      if(istatus.lt.1) then
         stop 'error returned from get_lapsdata_2d'
      endif

      source = comment(1:12)
      return
      end


