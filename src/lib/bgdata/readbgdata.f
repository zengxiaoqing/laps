      subroutine read_bgdata(nx_bg,ny_bg,
     +	   nzbg_ht,nzbg_tp,nzbg_sh,nzbg_uv,nzbg_ww,ctype
     +    ,bgpath,fname_bg,af_bg,fullname,cmodel,bgmodel
     +    ,prbght,prbgsh,prbguv,prbgww
     +    ,htbg, tpbg,uwbg,vwbg,shbg,wwbg
     +    ,htbg_sfc,prbg_sfc,shbg_sfc,tdbg_sfc,tpbg_sfc
     +    ,t_at_sfc,uwbg_sfc,vwbg_sfc,mslpbg,pcpbg,crefbg,tpw,cwat,swi
     +    ,istatus)

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
      integer i4_initial
      integer i4_valid
      integer i4time
      integer i4hr
      integer nan_flag
      integer nf_fid,nf_vid,nf_status
      integer istatus


c *** Background model grid data.
c
      real :: prbght(nx_bg,ny_bg,nzbg_ht) !Pressure (mb) ht and temp
      real :: prbgsh(nx_bg,ny_bg,nzbg_sh) !Pressure (mb) q
      real :: prbguv(nx_bg,ny_bg,nzbg_uv) !Pressure (mb) u- v-components
      real :: prbgww(nx_bg,ny_bg,nzbg_ww) !Pressure (mb) omega
      real  pr(nzbg_ht)

      real :: htbg(nx_bg,ny_bg,nzbg_ht)   !Height (m)
      real :: tpbg(nx_bg,ny_bg,nzbg_tp)   !Temperature (K)
      real :: shbg(nx_bg,ny_bg,nzbg_sh)   !Specific humidity (kg/kg)
      real :: uwbg(nx_bg,ny_bg,nzbg_uv)   !U-wind (m/s)
      real :: vwbg(nx_bg,ny_bg,nzbg_uv)   !V-wind (m/s)
      real :: wwbg(nx_bg,ny_bg,nzbg_ww)   !W-wind (pa/s)

      real :: mslpbg(nx_bg,ny_bg)         !mslp  (mb)
      real :: htbg_sfc(nx_bg,ny_bg)
      real :: prbg_sfc(nx_bg,ny_bg)
      real :: shbg_sfc(nx_bg,ny_bg)       !Specific humidity (kg/kg)
      real :: uwbg_sfc(nx_bg,ny_bg)
      real :: vwbg_sfc(nx_bg,ny_bg)
      real :: tdbg_sfc(nx_bg,ny_bg)
      real :: tpbg_sfc(nx_bg,ny_bg)
      real :: t_at_sfc(nx_bg,ny_bg)       !Skin/Ground Temperature
      real :: pcpbg(nx_bg,ny_bg)          !Precip at surface, ACPC (k/m^2)
      real :: crefbg(nx_bg,ny_bg)         !Composite Reflectivity
      real tpw(nx_bg,ny_bg)
      real cwat(nx_bg,ny_bg)
      real swi(nx_bg,ny_bg)

c     Local variables for the time being
      real r01(nx_bg,ny_bg)
      real llr(nx_bg,ny_bg)
      real s8a(nx_bg,ny_bg)

      real      lon0,lat1,lat2
      real      ssh2
      real      r_missing_data
     
      real      argmin,argmax

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

      write(6,*)' Subroutine read_bgdata: dims are ',nx_bg,ny_bg,nzbg_ht

      call get_r_missing_data(r_missing_data,istatus)

!     Initialize
      htbg_sfc = 0.
      prbg_sfc = r_missing_data
      tdbg_sfc = r_missing_data
      shbg_sfc = r_missing_data

      call s_len(cmodel,lencm)

      if(bgmodel.eq.0)then
        if(cmodel(1:lencm).eq.'LAPS_FUA'.or.
     +     cmodel(1:lencm).eq.'LAPS'.or.
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
     +                          ,nx_bg, ny_bg, nzbg_ht
     +                          ,htbg, pr, wwbg, shbg, tpbg, uwbg, vwbg       
     +                          ,uwbg_sfc, vwbg_sfc, tpbg_sfc, tdbg_sfc       
     +                          ,prbg_sfc, mslpbg, htbg_sfc, r01, pcpbg
     +                          ,crefbg, llr, s8a, swi, tpw
     +                          ,istatus)
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

         elseif(cmodel.eq.'LAPS_FUA')then

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
c
c ---------------Read LAPS Analyses---------------------
c
         elseif(cmodel.eq.'LAPS')then

           call get_laps_3d_analysis_data(i4_initial,nx_bg,ny_bg
     +,nzbg_ht, htbg,tpbg,uwbg,vwbg,shbg,wwbg,istatus)

           call get_laps_2d(i4_initial,'lsx','PS ',units_2d,comment_2d
     +,nx_bg,ny_bg,prbg_sfc,istatus)
           call get_laps_2d(i4_initial,'lsx','U  ',units_2d,comment_2d
     +,nx_bg,ny_bg,uwbg_sfc,istatus)
           call get_laps_2d(i4_initial,'lsx','V  ',units_2d,comment_2d
     +,nx_bg,ny_bg,vwbg_sfc,istatus)
           call get_laps_2d(i4_initial,'lsx','T  ',units_2d,comment_2d
     +,nx_bg,ny_bg,tpbg_sfc,istatus)
           call get_laps_2d(i4_initial,'lsx','TD ',units_2d,comment_2d
     +,nx_bg,ny_bg,shbg_sfc,istatus)
           call get_laps_2d(i4_initial,'lsx','MSL',units_2d,comment_2d
     +,nx_bg,ny_bg,mslpbg,istatus)

           if(istatus.ne.1)then
              print*,'Error returned: read_laps_analysis'
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
 
      elseif (bgmodel .eq. 5) then ! Process 40km RUC public-netcdf data

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
     .        bgmodel .eq. 8 .or.
     .        bgmodel .eq.12) then ! Process AVN, ECMWF or NOGAPS1.0 grib data

             call read_dgprep(bgmodel,cmodel,bgpath
     .,fname_bg,af_bg,nx_bg,ny_bg,nzbg_ht
     .,prbght,htbg,tpbg,shbg,uwbg,vwbg,wwbg
     .,htbg_sfc,prbg_sfc,shbg_sfc,tpbg_sfc,t_at_sfc
     .,uwbg_sfc,vwbg_sfc,mslpbg,istatus)

             if(istatus.ne.0)goto 99

             if(bgmodel.eq.12.and.
     &          cmodel(1:lencm).eq.'FMI_NETCDF_LL')then
                tdbg_sfc=shbg_sfc
             endif

             prbgsh=prbght
             prbguv=prbght
             prbgww=prbght

      elseif (bgmodel .eq. 9) then ! Process NWS Conus data (RUC,ETA,NGM,AVN)

             call read_conus_nws(bgpath,fname_bg,af_bg
     .,nx_bg,ny_bg,nzbg_ht,prbght
     .,htbg,tpbg,shbg,uwbg,vwbg
     .,gproj,lon0,lat1,lat2,istatus)
c
C WNI-BLS
      elseif (bgmodel .eq. 10) then ! Process Unidata NetCDF
       IF ( (cmodel .EQ. 'RUC_ISO') .OR.
     +      (cmodel .EQ. 'GFS_ISO') )THEN
          call read_unidata_iso(fullname,af_bg,cmodel
     .   ,nx_bg,ny_bg,nzbg_ht,nzbg_tp,nzbg_sh,nzbg_uv,nzbg_ww
     .   ,prbght,prbgsh,prbguv,prbgww
     .   ,htbg,tpbg,shbg,uwbg,vwbg,wwbg
     .   ,htbg_sfc,prbg_sfc,uwbg_sfc,vwbg_sfc,shbg_sfc
     .   ,tpbg_sfc,mslpbg,ctype,istatus)
       elseif(cmodel .EQ. 'RUC_HYB') THEN
           call read_unidata_ruc_hyb(fullname,af_bg,cmodel
     .   ,nx_bg,ny_bg,nzbg_ht,nzbg_tp,nzbg_sh,nzbg_uv,nzbg_ww
     .   ,prbght,prbgsh,prbguv,prbgww
     .   ,htbg,tpbg,shbg,uwbg,vwbg,wwbg
     .   ,htbg_sfc,prbg_sfc,uwbg_sfc,vwbg_sfc,shbg_sfc
     .   ,tpbg_sfc,mslpbg,ctype,istatus)
         print *, "Completed read of RUC_HYB"
       endif

      elseif (bgmodel .eq. 13) then ! Process GRIB1/GRIB2

         write(*,*) 'CALL DEGRIB_DATA: dims are ',nx_bg, ny_bg, nzbg_ht
         write(*,*) ' grib filename ',trim(fullname)

         call degrib_data(fullname, nx_bg, ny_bg, nzbg_ht, 
     &      prbght, htbg, tpbg, shbg, uwbg, vwbg, wwbg, 
     &      htbg_sfc, tpbg_sfc, shbg_sfc, uwbg_sfc, vwbg_sfc, 
     &      tdbg_sfc, t_at_sfc, prbg_sfc, mslpbg, pcpbg, crefbg, 
     &      tpw,cwat,istatus)

            prbgsh(:,:,:)=prbght(:,:,:) 
            prbguv(:,:,:)=prbght(:,:,:) 
            prbgww(:,:,:)=prbght(:,:,:) 

c        write(*, *) "READBGDATA htbg(3,30,1)", htbg(3,30,1)
c        write(*, *) "READBGDATA tpbg(3,30,1)", tpbg(3,30,1)
c        write(*, *) "READBGDATA wwbg(3,30,1)", wwbg(3,30,1)
c        write(*, *) "READBGDATA: shbg_sfc/tdbg_sfc is actually rh?"
         write(*, *) "READBGDATA shbg_sfc(3,30)", shbg_sfc(3,30)
         write(*, *) "READBGDATA tdbg_sfc(3,30)", tdbg_sfc(3,30)
c        do j = 1, nzbg_ht 
c           write(*, *) "READBGDATA pcpbg(3,",j, pcpbg(3,j)
c        enddo

      endif
c      
c - 3d fields
c
      call check_nan3 (htbg,nx_bg,ny_bg,nzbg_ht,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' ERROR: NaN found in height array '
c        goto 99
      endif

      call check_nan3 (tpbg,nx_bg,ny_bg,nzbg_tp,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' ERROR: NaN found in temperature array '
c        goto 99
      endif

      call check_nan3 (shbg,nx_bg,ny_bg,nzbg_sh,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' ERROR: NaN found in specif hum array '
c        goto 99
      endif

      call check_nan3 (uwbg,nx_bg,ny_bg,nzbg_uv,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' ERROR: NaN found in u-comp array '
c        goto 99
      endif

      call check_nan3 (vwbg,nx_bg,ny_bg,nzbg_uv,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' ERROR: NaN found in v-comp array '
c        goto 99
      endif

      call check_nan3 (wwbg,nx_bg,ny_bg,nzbg_ww,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' ERROR: NaN found in w-comp array '
      endif
c
c - 2d fields
c
      call check_nan2(htbg_sfc,nx_bg,ny_bg,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' ERROR: NaN found in sfc height array '
      endif

      call check_nan2(prbg_sfc,nx_bg,ny_bg,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' ERROR: NaN found in sfc pressure array '
      endif

      call check_nan2(tdbg_sfc,nx_bg,ny_bg,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' ERROR: NaN found in sfc temp array '
      endif

      call check_nan2(tpbg_sfc,nx_bg,ny_bg,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' ERROR: NaN found in sfc temp array '
      endif

      call check_nan2(shbg_sfc,nx_bg,ny_bg,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' ERROR: NaN found in sfc spec hum array '
      endif

      call check_nan2(uwbg_sfc,nx_bg,ny_bg,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' ERROR: NaN found in sfc u-comp array '
      endif

      call check_nan2(vwbg_sfc,nx_bg,ny_bg,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' ERROR: NaN found in sfc v-comp array '
      endif

      call check_nan2(mslpbg,nx_bg,ny_bg,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' ERROR: NaN found in sfc mslp array '
      endif

      call check_nan2(pcpbg,nx_bg,ny_bg,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' ERROR: NaN found in sfc pcpbg array '
      endif

      argmin = minval(prbg_sfc)
      argmax = maxval(prbg_sfc)
      write(6,*)' prbg_sfc range = ',argmin,argmax
      if(argmax .lt. 100. .or. argmax .gt. 1e6)then
          write(6,*)' WARNING: prbg_sfc has questionable range'
      endif

      write(6,*)' tdbg_sfc range = ',minval(tdbg_sfc),maxval(tdbg_sfc)
      write(6,*)' shbg_sfc range = ',minval(shbg_sfc),maxval(shbg_sfc)
      write(6,*)' Returning from read_bgdata'

      istatus = 0

99    if(istatus.ne. 0)then
         print*,'Error with background model data in read_bgdata'
      endif
      return
      end
