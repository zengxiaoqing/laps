      subroutine read_bgdata(nx_bg,ny_bg,
     +	   nzbg_ht,nzbg_tp,nzbg_sh,nzbg_uv,nzbg_ww,ctype
     +    ,bgpath,fname_bg,af_bg,fullname,cmodel,bgmodel
     +    ,prbght,prbgsh,prbguv,prbgww
     +    ,htbg, tpbg,uwbg,vwbg,shbg,wwbg
     +    ,htbg_sfc,prbg_sfc,shbg_sfc,tpbg_sfc
     +    ,uwbg_sfc,vwbg_sfc,mslpbg,istatus)

      implicit none
      include 'netcdf.inc'

      integer nx_bg
      integer ny_bg
      integer nzbg_ht
      integer nzbg_tp
      integer nzbg_sh
      integer nzbg_uv
      integer nzbg_ww
      integer bgmodel
      integer ntbg,nzbg
      integer i,j,k
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
      real, intent(out) :: tpbg_sfc(nx_bg,ny_bg)

      real      lon0,lat1,lat2

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

c     .   ,gproj,lon0,lat1,lat2,istatus)

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

      if(bgmodel.eq.0)then

        print*,'Reading fua/fsf'

        call get_directory_length(fullname,lend)
        directory=fullname(1:lend)
        call cv_asc_i4time(fullname(lend+1:lend+9),i4_initial)
        read(af_bg,100)i4hr
100     format(i4.4)

        i4_valid=i4_initial+i4hr/100*3600

        if(cmodel(1:lencm).eq.'LAPS_FUA')then

         call get_modelfg_3d(i4_valid,'U3 ',nx_bg,ny_bg,nzbg_uv
     .,uwbg,istatus)
         call get_modelfg_3d(i4_valid,'V3 ',nx_bg,ny_bg,nzbg_uv
     .,vwbg,istatus)
         call get_modelfg_3d(i4_valid,'T3 ',nx_bg,ny_bg,nzbg_ht
     .,tpbg,istatus)
         call get_modelfg_3d(i4_valid,'HT ',nx_bg,ny_bg,nzbg_ht
     .,htbg,istatus)
         call get_modelfg_3d(i4_valid,'SH ',nx_bg,ny_bg,nzbg_sh
     .,shbg,istatus)
         call get_modelfg_3d(i4_valid,'OM ',nx_bg,ny_bg,nzbg_ww
     .,wwbg,istatus)

c     print*,'get pressures from pressure.nl'
         call get_pres_1d(i4time,nzbg_ht,pr,istatus)
         do k = 1,nzbg_ht
            pr(k)=pr(k)/100.
            do j=1,ny_bg
            do i=1,nx_bg
               prbght(i,j,k)=pr(k)
               prbgsh(i,j,k)=pr(k)
               prbguv(i,j,k)=pr(k)
               prbgww(i,j,k)=pr(k)
            enddo
            enddo
         enddo

         call get_modelfg_2d(i4_valid,'USF',nx_bg,ny_bg,uwbg_sfc
     .,istatus)
         call get_modelfg_2d(i4_valid,'VSF',nx_bg,ny_bg,vwbg_sfc
     .,istatus)
         call get_modelfg_2d(i4_valid,'TSF',nx_bg,ny_bg,tpbg_sfc
     .,istatus)
         call get_modelfg_2d(i4_valid,'DSF',nx_bg,ny_bg,shbg_sfc
     .,istatus)
         call get_modelfg_2d(i4_valid,'SLP',nx_bg,ny_bg,mslpbg
     .,istatus)
         call get_modelfg_2d(i4_valid,'PSF',nx_bg,ny_bg,prbg_sfc
     .,istatus)

        elseif(cmodel.eq.'MODEL_FUA')then

         nf_status = NF_OPEN(fullname,NF_NOWRITE,nf_fid)
         if(nf_status.ne.NF_NOERR) then
            print *, NF_STRERROR(nf_status)
            print *,'NF_OPEN ', fullname
            return
         endif
         nf_status=NF_INQ_VARID(nf_fid,'z',nf_vid)
         if(nf_status.ne.NF_NOERR) then
            print *, NF_STRERROR(nf_status)
            print *,'in NF_GET_VAR_ model '
            return
         endif
         nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,nzbg)
         if(nf_status.ne.NF_NOERR) then
            print *, NF_STRERROR(nf_status)
            print *,'dim n_valtimes'
            return
         endif
         nf_status=NF_INQ_VARID(nf_fid,'levels',nf_vid)
         if(nf_status.ne.NF_NOERR) then
            print *, NF_STRERROR(nf_status)
            print *,'in NF_GET_VAR_ model '
            return
         endif
         nf_status=NF_GET_VARA_INT(nf_fid,nf_vid,1,nzbg_ht,pr)
         if(nf_status.ne.NF_NOERR) then
            print *, NF_STRERROR(nf_status)
            print *,'in NF_GET_VAR_ model '
            return
         endif

         do k = 1,nzbg_ht
            do j=1,ny_bg
            do i=1,nx_bg
               prbght(i,j,k)=pr(k)
               prbgsh(i,j,k)=pr(k)
               prbguv(i,j,k)=pr(k)
               prbgww(i,j,k)=pr(k)
            enddo
            enddo
         enddo

c upper air
         call get_lapsdata_3d(i4_initial,i4_valid,nx_bg
     1            ,ny_bg,nzbg_uv,directory,'U3 '
     1            ,units_2d,comment_2d,uwbg,istatus)

         call get_lapsdata_3d(i4_initial,i4_valid,nx_bg
     1            ,ny_bg,nzbg_uv,directory,'V3 '
     1            ,units_2d,comment_2d,vwbg,istatus)

         call get_lapsdata_3d(i4_initial,i4_valid,nx_bg
     1            ,ny_bg,nzbg_ht,directory,'T3 '
     1            ,units_2d,comment_2d,tpbg,istatus)

         call get_lapsdata_3d(i4_initial,i4_valid,nx_bg
     1            ,ny_bg,nzbg_ht,directory,'HT '
     1            ,units_2d,comment_2d,htbg,istatus)

         call get_lapsdata_3d(i4_initial,i4_valid,nx_bg
     1            ,ny_bg,nzbg_sh,directory,'SH '
     1            ,units_2d,comment_2d,shbg,istatus)

         call get_lapsdata_3d(i4_initial,i4_valid,nx_bg
     1            ,ny_bg,nzbg_ww,directory,'OM '
     1            ,units_2d,comment_2d,wwbg,istatus)

c sfc data
         call get_lapsdata_2d(i4time,i4_valid,directory,'SLP'
     1            ,units_2d,comment_2d,nx_bg,ny_bg,mslpbg
     1            ,istatus)

         call get_lapsdata_2d(i4time,i4_valid,directory,'DSF'
     1            ,units_2d,comment_2d,nx_bg,ny_bg,shbg_sfc
     1            ,istatus)

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
     .,htbg_sfc,prbg_sfc,shbg_sfc,tpbg_sfc
     .,uwbg_sfc,vwbg_sfc,mslpbg,istatus)

          if(istatus.ne.0)goto 99

          if(ctype.eq."lapsb")then

c convert rh to sh.
             print*,'Prepare grids for laps background-lga'
             call lprep_eta_conusc(nx_bg,ny_bg,nzbg_ht
     +,prbght,tpbg,shbg,tpbg_sfc,prbg_sfc,shbg_sfc,istatus)
 
          endif

          prbgsh=prbght
          prbguv=prbght
          prbgww=prbght
c
      elseif (bgmodel .eq. 4) then ! Process SBN Conus 211 data (Eta or RUC)

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
