      subroutine read_bgdata(nx_bg,ny_bg,
     +	   nzbg_ht,nzbg_tp,nzbg_sh,nzbg_uv,nzbg_ww,ctype
     +    ,bgpath,fname_bg,af_bg,fullname,cmodel,bgmodel
     +    ,prbght,prbgsh,prbguv,prbgww
     +    ,htbg, tpbg,uwbg,vwbg,shbg,wwbg
     +    ,htbg_sfc,prbg_sfc,shbg_sfc,tdbg_sfc,tpbg_sfc
     +    ,uwbg_sfc,vwbg_sfc,mslpbg,istatus)

c KML: CHANGES MADE APRIL 2004
c tdbg_sfc (model 2m dew point) is now read in during subroutine read_eta_conusc
c tdbg_sfc array is checked for nan
c KML: END

      implicit none
      include 'netcdf.inc'

      integer nx_bg
      integer ny_bg
      integer nx_l
      integer ny_l
      integer nzbg_ht
      integer nzbg_tp
      integer nzbg_sh
      integer nzbg_uv
      integer nzbg_ww
      integer bgmodel
      integer ntbg,nzbg
      integer i,j,k,l
      integer lencm,lend
      integer ivaltimes(100)
      integer i4_initial
      integer i4_valid
      integer i4time
      integer i4hr
      integer nan_flag
      integer nf_fid,nf_vid,nf_status
      integer istatus
c
c *** Background model grid data.
c
      real, intent(out) :: prbght(nx_bg,ny_bg,nzbg_ht) !Pressure (mb) ht and temp
      real, intent(out) :: prbgsh(nx_bg,ny_bg,nzbg_sh) !Pressure (mb) q
      real, intent(out) :: prbguv(nx_bg,ny_bg,nzbg_uv) !Pressure (mb) u- v-components
      real, intent(out) :: prbgww(nx_bg,ny_bg,nzbg_ww) !Pressure (mb) omega
      real  pr(nzbg_ht)

      real, intent(out) :: htbg(nx_bg,ny_bg,nzbg_ht)   !Height (m)
      real, intent(out) :: tpbg(nx_bg,ny_bg,nzbg_tp)   !Temperature (K)
      real, intent(out) :: shbg(nx_bg,ny_bg,nzbg_sh)   !Specific humidity (kg/kg)
      real, intent(out) :: uwbg(nx_bg,ny_bg,nzbg_uv)   !U-wind (m/s)
      real, intent(out) :: vwbg(nx_bg,ny_bg,nzbg_uv)   !V-wind (m/s)
      real, intent(out) :: wwbg(nx_bg,ny_bg,nzbg_ww)   !W-wind (pa/s)

      real, intent(out) :: mslpbg(nx_bg,ny_bg)         !mslp  (mb)
      real, intent(out) :: htbg_sfc(nx_bg,ny_bg)
      real, intent(out) :: prbg_sfc(nx_bg,ny_bg)
      real, intent(out) :: shbg_sfc(nx_bg,ny_bg)
      real, intent(out) :: uwbg_sfc(nx_bg,ny_bg)
      real, intent(out) :: vwbg_sfc(nx_bg,ny_bg)
      real, intent(out) :: tdbg_sfc(nx_bg,ny_bg)
      real, intent(out) :: tpbg_sfc(nx_bg,ny_bg)

      real      lon0,lat1,lat2
      real      ssh2
      real      r_missing_data

      character*256   bgpath
      character*256   fname_bg
      character*200   fullname
      character*150   directory
      character*132   cmodel
      character*125   comment_2d
      character*10    units_2d
      character*13    fname13
      character*13    fname9_to_wfo_fname13
      character*6     c6_maproj
      character*5     ctype      !either "dprep" or "lapsb" depending on dprep or lga
      character*4     af_bg
      character*2     gproj


      interface

         subroutine read_sbn_grids(cdfname,af,cmodel,
     .nxbg,nybg,nzbght,nzbgtp,nzbgsh,nzbguv,nzbgww,
     .prbght,prbgsh,prbguv,prbgww,
     .ht,tp,sh,uw,vw,ww,
     .ht_sfc,pr_sfc,uw_sfc,vw_sfc,sh_sfc,tp_sfc,mslp,
     .ctype,istatus)

         real  ::   pr_sfc(nxbg,nybg)
         real  ::   uw_sfc(nxbg,nybg)
         real  ::   vw_sfc(nxbg,nybg)
         real  ::   sh_sfc(nxbg,nybg)
         real  ::   tp_sfc(nxbg,nybg)
         real  ::   ht_sfc(nxbg,nybg)
         real  ::     mslp(nxbg,nybg)
c
         real  :: prbght(nxbg,nybg,nzbght)
         real  :: prbgsh(nxbg,nybg,nzbgsh)
         real  :: prbguv(nxbg,nybg,nzbguv)
         real  :: prbgww(nxbg,nybg,nzbgww)
         real  ::     ht(nxbg,nybg,nzbght)
         real  ::     tp(nxbg,nybg,nzbgtp)
         real  ::     sh(nxbg,nybg,nzbgsh)
         real  ::     uw(nxbg,nybg,nzbguv)
         real  ::     vw(nxbg,nybg,nzbguv)
         real  ::     ww(nxbg,nybg,nzbgww)

         character*200 cdfname
         character*132 cmodel
         character*5   ctype
         character*4   af
         integer       nxbg
         integer       nybg
         integer       nzbght
         integer       nzbgsh
         integer       nzbguv
         integer       nzbgww
         integer       istatus
c
         end subroutine

         subroutine read_dgprep(bgmodel,cmodel,path,fname,af
     .   ,nx,ny,nz
     .   ,pr,ht,tp,sh,uw,vw,ww
     .   ,ht_sfc,pr_sfc,td_sfc,tp_sfc
     .   ,uw_sfc,vw_sfc,mslp,istatus)

         integer bgmodel
         integer nx,ny,nz
         integer istatus

         real*4 ht(nx,ny,nz)
         real*4 tp(nx,ny,nz)
         real*4 sh(nx,ny,nz)
         real*4 uw(nx,ny,nz)
         real*4 vw(nx,ny,nz)
         real*4 ww(nx,ny,nz)
         real*4 pr(nx,ny,nz)

         real*4 ht_sfc(nx,ny)
         real*4 pr_sfc(nx,ny)
         real*4 td_sfc(nx,ny)
         real*4 tp_sfc(nx,ny)
         real*4 uw_sfc(nx,ny)
         real*4 vw_sfc(nx,ny)
         real*4 mslp(nx,ny)

c        real*4 lon0
c        real*4 lat1,lat2

         character cmodel*132
         character path*256
         character fname*200
         character af*4
c        character gproj*2
 
         end subroutine

      end interface

      call s_len(cmodel,lencm)

      call get_r_missing_data(r_missing_data,istatus)

      if(bgmodel.eq.0)then
        if(cmodel(1:lencm).eq.'LAPS_FUA'.or.
     +     cmodel(1:lencm).eq.'MODEL_FUA')then

           call get_directory_length(fullname,lend)
           directory=fullname(1:lend)
           call cv_asc_i4time(fullname(lend+1:lend+9),i4_initial)
           read(af_bg,100)i4hr
100     format(i4.4)
           i4_valid=i4_initial+i4hr/100*3600

c the following subroutine should also work for different
c domain fua/fsf but we'll try the get_lapsdata stuff first.

         if(cmodel(1:lencm).eq.'MODEL_FUA')then
            call read_fuafsf_cdf(fullname
     +,nx_bg, ny_bg, nzbg_ht
     +,htbg, pr, wwbg, shbg, tpbg, uwbg, vwbg
     +,uwbg_sfc, vwbg_sfc, tpbg_sfc, shbg_sfc
     +,prbg_sfc, mslpbg, htbg_sfc, istatus)
            if(istatus.ne.1)then
               print*,'Error returned: read_fuafsf_cdf'
               return
            endif
            do j=1,ny_bg
            do i=1,nx_bg
               prbght(i,j,:)=pr(:)
               prbgsh(i,j,:)=pr(:)
               prbguv(i,j,:)=pr(:)
               prbgww(i,j,:)=pr(:)
            enddo
            enddo

         else  ! cmodel = LAPS_FUA: Same domain!

           call get_grid_dim_xy(nx_l,ny_l,istatus)

           if(nx_l .ne. nx_bg .or. ny_l .ne. ny_bg)then
              print*
              print*,' ***********************************'
              print*,' Error: nx-laps ne nx_bg '
              print*,' background.nl var cmodel = ',cmodel
              print*,' ***********************************'
              print*
              return
           endif

           call s_len(fullname,lend)
           fullname=fullname(1:lend)//".fua"
           nf_status = NF_OPEN(fullname,NF_NOWRITE,nf_fid)
           if(nf_status.ne.NF_NOERR) then
              print *, NF_STRERROR(nf_status)
              print *,'NF_OPEN ', fullname
              return
           endif

c        nf_status=NF_INQ_VARID(nf_fid,'z',nf_vid)
c        if(nf_status.ne.NF_NOERR) then
c           print *, NF_STRERROR(nf_status)
c           print *,'in NF_GET_VAR_ model '
c           return
c        endif
c        nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,nzbg)
c        if(nf_status.ne.NF_NOERR) then
c           print *, NF_STRERROR(nf_status)
c           print *,'dim n_valtimes'
c           return
c        endif

           nf_status=NF_INQ_VARID(nf_fid,'level',nf_vid)
           if(nf_status.ne.NF_NOERR) then
              print *, NF_STRERROR(nf_status)
              print *,'in NF_GET_VAR_ model '
              return
           endif
           nf_status=NF_GET_VAR_REAL(nf_fid,nf_vid,pr)
           if(nf_status.ne.NF_NOERR) then
              print *, NF_STRERROR(nf_status)
              print *,'in NF_GET_VAR_ model '
              return
           endif

           do k = 1,nzbg_ht
           do j=1,ny_bg
           do i=1,nx_bg
              prbght(i,j,k)=pr(nzbg_ht-k+1)
              prbgsh(i,j,k)=pr(nzbg_ht-k+1)
              prbguv(i,j,k)=pr(nzbg_ht-k+1)
              prbgww(i,j,k)=pr(nzbg_ht-k+1)
           enddo
           enddo
           enddo

c upper air
           call get_lapsdata_3d(i4_initial,i4_valid,nx_bg
     1            ,ny_bg,nzbg_uv,directory,'U3 '
     1            ,units_2d,comment_2d,uwbg,istatus)
           if(istatus.ne.1)then
              print*,'Error 3D bkgd file (U3): ',directory(1:lend)
              return
           endif

           call get_lapsdata_3d(i4_initial,i4_valid,nx_bg
     1            ,ny_bg,nzbg_uv,directory,'V3 '
     1            ,units_2d,comment_2d,vwbg,istatus)
           if(istatus.ne.1)then
              print*,'Error 3D bkgd file (V3): ',directory(1:lend)
              return
           endif

           call get_lapsdata_3d(i4_initial,i4_valid,nx_bg
     1            ,ny_bg,nzbg_ht,directory,'T3 '
     1            ,units_2d,comment_2d,tpbg,istatus)
           if(istatus.ne.1)then
              print*,'Error 3D bkgd file (T3): ',directory(1:lend)
              return
           endif

           call get_lapsdata_3d(i4_initial,i4_valid,nx_bg
     1            ,ny_bg,nzbg_ht,directory,'HT '
     1            ,units_2d,comment_2d,htbg,istatus)
           if(istatus.ne.1)then
              print*,'Error 3D bkgd file (HT): ',directory(1:lend)
              return
           endif

           call get_lapsdata_3d(i4_initial,i4_valid,nx_bg
     1            ,ny_bg,nzbg_sh,directory,'SH '
     1            ,units_2d,comment_2d,shbg,istatus)
           if(istatus.ne.1)then
              print*,'Error 3D bkgd file (SH): ',directory(1:lend)
              return
           endif

           call get_lapsdata_3d(i4_initial,i4_valid,nx_bg
     1            ,ny_bg,nzbg_ww,directory,'OM '
     1            ,units_2d,comment_2d,wwbg,istatus)

           if(istatus.ne.1)then
              print*,'Error 3D bkgd file (OM): ',directory(1:lend)
              return
           endif
c sfc data
           search_dir: do l=lend,1,-1
              if(directory(l:l).eq.'f')then
                 if(directory(l:l+2).eq.'fua')then
                    exit search_dir
                 endif
              endif
           enddo search_dir

           if(l.le.1)then
              print*,'Unable to determine location of fua in'
              print*,'directory string for fsf. Return with no data'
              return
           endif

           directory(l:l+2)='fsf'

           call get_lapsdata_2d(i4_initial,i4_valid,directory
     1            ,'PSF',units_2d,comment_2d,nx_bg,ny_bg,prbg_sfc
     1            ,istatus)
           if(istatus.ne.1)then
              print*,'Error 2D bkgd file (PSF): ',directory(1:lend)
              return
           endif

c          prbg_sfc=prbg_sfc*100.

           call get_lapsdata_2d(i4_initial,i4_valid,directory
     1            ,'TSF',units_2d,comment_2d,nx_bg,ny_bg,tpbg_sfc
     1            ,istatus)
           if(istatus.ne.1)then
              print*,'Error 2D bkgd file (TSF): ',directory(1:lend)
              return
           endif

           call get_lapsdata_2d(i4_initial,i4_valid,directory
     1            ,'SLP',units_2d,comment_2d,nx_bg,ny_bg,mslpbg
     1            ,istatus)
           if(istatus.ne.1)then
              print*,'Error 2D bkgd file (SLP): ',directory(1:lend)
              return
           endif

           call get_lapsdata_2d(i4_initial,i4_valid,directory
     1            ,'DSF',units_2d,comment_2d,nx_bg,ny_bg,shbg_sfc
     1            ,istatus)
           if(istatus.ne.1)then
              print*,'Error 2D bkgd file (DSF): ',directory(1:lend)
              return
           endif

           call get_lapsdata_2d(i4_initial,i4_valid,directory
     1            ,'USF',units_2d,comment_2d,nx_bg,ny_bg,uwbg_sfc
     1            ,istatus)
           if(istatus.ne.1)then
              print*,'Error 2D bkgd file (USF): ',directory(1:lend)
              return
           endif

           call get_lapsdata_2d(i4_initial,i4_valid,directory
     1            ,'VSF',units_2d,comment_2d,nx_bg,ny_bg,vwbg_sfc
     1            ,istatus)
           if(istatus.ne.1)then
              print*,'Error 2D bkgd file (VSF): ',directory(1:lend)
              return
           endif

           call get_lapsdata_2d(i4_initial,i4_valid,directory
     1            ,'TER',units_2d,comment_2d,nx_bg,ny_bg,htbg_sfc
     1            ,istatus)
           if(istatus.ne.1)then
              print*,'Error 2D bkgd file (TER): ',directory(1:lend)
              return
           endif

         endif !MODEL_FUA?!

        endif

      elseif (bgmodel .eq. 1) then     ! Process 60 km RUC data

          call read_ruc60_native(bgpath,fname_bg,af_bg,nx_bg,ny_bg
     .,nzbg_ht,prbght,htbg,tpbg,shbg,uwbg,vwbg,gproj,lon0,istatus)

      elseif (bgmodel .eq. 2) then ! Process 48 km ETA conus-c grid data
c
c for now all fields have 3D dimension of nzbg_ht
c
          call read_eta_conusC(fullname,nx_bg,ny_bg,nzbg_ht
     .,htbg,prbght,tpbg,uwbg,vwbg,shbg,wwbg
     .,htbg_sfc,prbg_sfc,shbg_sfc,tdbg_sfc,tpbg_sfc
     .,uwbg_sfc,vwbg_sfc,mslpbg,istatus)

          if(istatus.ne.0)goto 99

          if(ctype.eq."lapsb")then

c convert rh to sh.
             print*,'Prepare grids for laps background-lga'
             call lprep_eta_conusc(nx_bg,ny_bg,nzbg_ht
     +,prbght,tpbg,shbg,tpbg_sfc,prbg_sfc,shbg_sfc,istatus)
 
          endif
          if(cmodel(1:lencm).eq.'ORSM_HKO')then
             print*,'In read_bgdata'
             print*,'Convert ww to pa/sec for ORSM_HKO'
             wwbg=wwbg/36.
          endif

          prbgsh=prbght
          prbguv=prbght
          prbgww=prbght
c
      elseif (bgmodel .eq. 4) then

          call read_sbn_grids(fullname,af_bg,cmodel
     .,nx_bg,ny_bg,nzbg_ht,nzbg_tp,nzbg_sh,nzbg_uv,nzbg_ww
     .,prbght,prbgsh,prbguv,prbgww
     .,htbg,tpbg,shbg,uwbg,vwbg,wwbg
     .,htbg_sfc,prbg_sfc,uwbg_sfc,vwbg_sfc,shbg_sfc
     .,tpbg_sfc,mslpbg,ctype,istatus)

          if(istatus.ne.0)goto 99

c original nz allows reading of sfc variables (ie., 2m and 10m)
c but now decrement nz since arrays only contain 3d info.

c         if(cmodel(1:lencm).eq.'RUC40_NATIVE')then
c            nzbg_sh=nzbg_sh-1
c            nzbg_uv=nzbg_uv-1
c         elseif(cmodel(1:lencm).eq.'ETA48_CONUS')then
c            nzbg_sh=nzbg_sh-1
c            nzbg_uv=nzbg_uv-1
c         elseif(cmodel(1:lencm).eq.'AVN_SBN_CYLEQ')then
c            nzbg_sh=nzbg_sh-2
c            nzbg_uv=nzbg_uv-2
c         endif
 
      elseif (bgmodel .eq. 5) then ! Process 40 km RUC data

          call read_ruc2_hybb(fullname,nx_bg,ny_bg,nzbg_ht
     +,mslpbg,htbg,prbght,shbg,uwbg,vwbg,tpbg,wwbg,istatus)

          if(istatus.ne.0)goto 99 

          print*,'Read complete'

          if(ctype.eq.'lapsb')then

             print*,' entering lprep_ruc2_hybrid'

             call lprep_ruc2_hybrid(nx_bg,ny_bg,nzbg_ht
     +,htbg,prbght,shbg,uwbg,vwbg,tpbg,uwbg_sfc,vwbg_sfc
     +,tpbg_sfc,prbg_sfc,shbg_sfc,htbg_sfc,istatus)
             print*,'Data prep complete'

          endif

          prbgsh=prbght
          prbguv=prbght
          prbgww=prbght
c
c ETA grib ingest currently disabled (J. Smart 9-4-98)
c Also, NOGAPS 2.5 degree obsolete.
c bgmodel 3 = FA (Taiwan). bgmodel 6 = NOGAPS1.0. bgmodel 8 = AVN 1.0 deg
c
      elseif (bgmodel .eq. 3 .or.
     .        bgmodel .eq. 6 .or.
     .        bgmodel .eq. 8) then ! Process AVN or NOGAPS1.0 grib data

             call read_dgprep(bgmodel,cmodel,bgpath
     .,fname_bg,af_bg,nx_bg,ny_bg,nzbg_ht
     .,prbght,htbg,tpbg,shbg,uwbg,vwbg,wwbg
     .,htbg_sfc,prbg_sfc,shbg_sfc,tpbg_sfc
     .,uwbg_sfc,vwbg_sfc,mslpbg,istatus)

             if(istatus.ne.0)goto 99

             prbgsh=prbght
             prbguv=prbght
             prbgww=prbght

      elseif (bgmodel .eq. 9) then ! Process NWS Conus data (RUC,ETA,NGM,AVN)

             call read_conus_nws(bgpath,fname_bg,af_bg
     .,nx_bg,ny_bg,nzbg_ht,prbght
     .,htbg,tpbg,shbg,uwbg,vwbg
     .,gproj,lon0,lat1,lat2,istatus)
c
      endif
      
      call checknan_3d(htbg,nx_bg,ny_bg,nzbg_ht,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' ERROR: NaN found in height array '
         goto 99
      endif

      call checknan_3d(tpbg,nx_bg,ny_bg,nzbg_tp,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' ERROR: NaN found in temperature array '
         goto 99
      endif

      call checknan_3d(shbg,nx_bg,ny_bg,nzbg_sh,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' ERROR: NaN found in specif hum array '
         goto 99
      endif

      call checknan_3d(uwbg,nx_bg,ny_bg,nzbg_uv,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' ERROR: NaN found in u-comp array '
         goto 99
      endif

      call checknan_3d(vwbg,nx_bg,ny_bg,nzbg_uv,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' ERROR: NaN found in v-comp array '
         goto 99
      endif

      call checknan_3d(wwbg,nx_bg,ny_bg,nzbg_ww,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' ERROR: NaN found in w-comp array '
      endif

      call checknan_2d(htbg_sfc,nx_bg,ny_bg,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' ERROR: NaN found in sfc height array '
      endif

      call checknan_2d(prbg_sfc,nx_bg,ny_bg,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' ERROR: NaN found in sfc pressure array '
      endif

      call checknan_2d(tdbg_sfc,nx_bg,ny_bg,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' ERROR: NaN found in sfc temp array '
      endif

      call checknan_2d(tpbg_sfc,nx_bg,ny_bg,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' ERROR: NaN found in sfc temp array '
      endif

      call checknan_2d(shbg_sfc,nx_bg,ny_bg,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' ERROR: NaN found in sfc spec hum array '
      endif

      call checknan_2d(uwbg_sfc,nx_bg,ny_bg,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' ERROR: NaN found in sfc u-comp array '
      endif

      call checknan_2d(vwbg_sfc,nx_bg,ny_bg,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' ERROR: NaN found in sfc v-comp array '
      endif

      call checknan_2d(mslpbg,nx_bg,ny_bg,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' ERROR: NaN found in sfc mslp array '
      endif
      istatus = 0

99    if(istatus.ne. 0)then
         print*,'Error with background model data in read_bgdata'
      endif
      return
      end
