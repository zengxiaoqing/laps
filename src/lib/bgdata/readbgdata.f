      subroutine read_bgdata(mxlvls,nx_bg,ny_bg,
     +	   nzbg_ht,nzbg_sh,nzbg_uv,nzbg_ww,ctype
     +    ,bgpath,fname_bg,af_bg,fullname,cmodel,bgmodel
     +    ,prbght,prbgsh,prbguv,prbgww
     +    ,htbg, tpbg,uwbg,vwbg,shbg,wwbg
     +    ,htbg_sfc,prbg_sfc,shbg_sfc,tpbg_sfc
     +    ,uwbg_sfc,vwbg_sfc,mslpbg,istatus)

      implicit none

      integer nx_bg
      integer ny_bg
      integer nzbg_ht
      integer nzbg_sh
      integer nzbg_uv
      integer nzbg_ww
      integer mxlvls
      integer bgmodel
      integer ntbg
      integer k
      integer lencm
      integer ivaltimes(100)
      integer istatus

c
c *** Background model grid data.
c
      real, intent(out) :: prbght(:,:,:)     !Pressure (mb) ht and temp
      real, intent(out) :: prbgsh(:,:,:)     !Pressure (mb) q
      real, intent(out) :: prbguv(:,:,:)     !Pressure (mb) u- v-components
      real, intent(out) :: prbgww(:,:,:)     !Pressure (mb) omega

      real, intent(out) :: htbg(:,:,:)     !Height (m)
      real, intent(out) :: tpbg(:,:,:)     !Temperature (K)
      real, intent(out) :: shbg(:,:,:)     !Specific humidity (kg/kg)
      real, intent(out) :: uwbg(:,:,:)     !U-wind (m/s)
      real, intent(out) :: vwbg(:,:,:)     !V-wind (m/s)
      real, intent(out) :: wwbg(:,:,:)     !W-wind (pa/s)
      real, intent(out) :: mslpbg(:,:)     !mslp  (mb)
      real, intent(out) :: htbg_sfc(:,:)
      real, intent(out) :: prbg_sfc(:,:)
      real, intent(out) :: shbg_sfc(:,:)
      real, intent(out) :: uwbg_sfc(:,:)
      real, intent(out) :: vwbg_sfc(:,:)
      real, intent(out) :: tpbg_sfc(:,:)

      real      lon0,lat1,lat2

      character*5     ctype      !either "dprep" or "lapsb" depending on dprep or lga
      character*2     gproj
      character*132   cmodel
      character*256   fullname
      character*256   bgpath
      character*256   fname_bg
      character*4     af_bg
      character*6     c6_maproj
      character*13    fname13
      character*13    fname9_to_wfo_fname13

      interface
         subroutine read_sbn_grids(cdfname,af,cmodel,
     .mxlvls,nxbg,nybg,nzbght,nzbgsh,nzbguv,nzbgww,
     .prbght,prbgsh,prbguv,prbgww,
     .ht,tp,sh,uw,vw,ww,
     .ht_sfc,pr_sfc,uw_sfc,vw_sfc,sh_sfc,tp_sfc,mslp,
     .ctype,istatus)

         real  ::   pr_sfc(:,:)
         real  ::   uw_sfc(:,:)
         real  ::   vw_sfc(:,:)
         real  ::   sh_sfc(:,:)
         real  ::   tp_sfc(:,:)
         real  ::   ht_sfc(:,:)
         real  ::     mslp(:,:)
c
         real  :: prbght(:,:,:)
         real  :: prbgsh(:,:,:)
         real  :: prbguv(:,:,:)
         real  :: prbgww(:,:,:)
         real  ::     ht(:,:,:)
         real  ::     tp(:,:,:)
         real  ::     sh(:,:,:)
         real  ::     uw(:,:,:)
         real  ::     vw(:,:,:)
         real  ::     ww(:,:,:)

         character*256 cdfname
         character*132 cmodel
         character*5   ctype
         character*4   af
         integer       mxlvls
         integer       nxbg
         integer       nybg
         integer       nzbght
         integer       nzbgsh
         integer       nzbguv
         integer       nzbgww
         integer       istatus
c
         end subroutine
      end interface

      call s_len(cmodel,lencm)

      if(bgmodel.eq.0)then

         print*,'readbgdata not yet designed to read laps -'
         print*,'this is temporary and will be added soon'

      elseif (bgmodel .eq. 1) then     ! Process 60 km RUC data
          call read_ruc60_native(bgpath,fname_bg,af_bg,nx_bg,ny_bg,
     .               nzbg_ht,nzbg_sh,nzbg_uv,nzbg_ww,
     .               prbght,htbg,tpbg,shbg,uwbg,vwbg,gproj,lon0,istatus)

      elseif (bgmodel .eq. 2) then ! Process 48 km ETA conus-c grid data
c
c for now all fields have 3D dimension of nzbg_ht
c
          call read_eta_conusC(fullname,nx_bg,ny_bg,nzbg_ht,
     .               htbg,prbght,tpbg,uwbg,vwbg,shbg,wwbg,
     .               htbg_sfc,prbg_sfc,shbg_sfc,tpbg_sfc,
     .               uwbg_sfc,vwbg_sfc,mslpbg,
     .               istatus)

          if(istatus.ne.0)goto 99

          if(ctype.eq."lapsb")then

c convert rh to sh.
             print*,'Preparing grids for laps background'
             call lprep_eta_conusc(nx_bg,ny_bg,nzbg_ht
     +,prbght,tpbg,shbg,tpbg_sfc,prbg_sfc,shbg_sfc,istatus)
 
          endif

          prbgsh=prbght
          prbguv=prbght
          prbgww=prbght
c
      elseif (bgmodel .eq. 4) then ! Process SBN Conus 211 data (Eta or RUC)

          print*,'entering read_sbn_grids'
          call read_sbn_grids(fullname,af_bg,cmodel,mxlvls,
     .         nx_bg,ny_bg,nzbg_ht,nzbg_sh,nzbg_uv,nzbg_ww,
     .         prbght,prbgsh,prbguv,prbgww,
     .         htbg,tpbg,shbg,uwbg,vwbg,wwbg,
     .         htbg_sfc,prbg_sfc,uwbg_sfc,vwbg_sfc,shbg_sfc,
     .         tpbg_sfc,mslpbg,ctype,istatus)

          if(istatus.ne.0)goto 99

c original nz allows reading of sfc variables (ie., 2m and 10m)
c but now decrement nz since arrays only contain 3d info.

          if(cmodel(1:lencm).eq.'RUC40_NATIVE')then
             nzbg_sh=nzbg_sh-1
             nzbg_uv=nzbg_uv-1
          elseif(cmodel(1:lencm).eq.'ETA48_CONUS')then
             nzbg_sh=nzbg_sh-1
             nzbg_uv=nzbg_uv-1
          elseif(cmodel(1:lencm).eq.'AVN_SBN_CYLEQ')then
             nzbg_sh=nzbg_sh-2
             nzbg_uv=nzbg_uv-2
          endif
 
      elseif (bgmodel .eq. 5) then ! Process 40 km RUC data

          call read_ruc2_hybb(fullname,nx_bg,ny_bg,nzbg_ht,
     +         mslpbg,htbg,prbght,shbg,uwbg,vwbg,tpbg,wwbg,istatus)

          if(istatus.ne.0)goto 99 

          print*,'Read complete'

          if(ctype.eq.'lapsb')then

             print*,' entering prep'
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
     .                 ,fname_bg,af_bg,nx_bg,ny_bg,nzbg_ht
     .                 ,prbght,htbg,tpbg,shbg,uwbg,vwbg,wwbg
     .                 ,htbg_sfc,prbg_sfc,shbg_sfc,tpbg_sfc
     .                 ,uwbg_sfc,vwbg_sfc,mslpbg
     .                 ,gproj,lon0,lat1,lat2,istatus)

             if(istatus.ne.0)goto 99

             prbgsh=prbght
             prbguv=prbght
             prbgww=prbght

      elseif (bgmodel .eq. 9) then ! Process NWS Conus data (RUC,ETA,NGM,AVN)
             call read_conus_nws(bgpath,fname_bg,af_bg,
     .               nx_bg,ny_bg,nzbg_ht,prbght
     .              ,htbg,tpbg,shbg,uwbg,vwbg,
     .               gproj,lon0,lat1,lat2,istatus)
c
      endif

99    if(istatus.ne. 0)then
         print*,'Error with background model data in read_bgdata'
      endif

      return
      end
