      subroutine read_bgdata(nx_bg,ny_bg,nz_bg
     +    ,bgpath,fname_bg,af_bg,fullname,cmodel,bgmodel
     +    ,htbg, prbg,tpbg,uwbg,vwbg,shbg,wwbg
     +    ,htbg_sfc,prbg_sfc,shbg_sfc,tpbg_sfc
     +    ,uwbg_sfc,vwbg_sfc,mslpbg
     +    ,gproj,lon0,lat1,lat2,istatus)

      implicit none

      integer nx_bg,ny_bg,nz_bg
      integer bgmodel
      integer nxbg,nybg,nzbg(5),ntbg
      integer istatus

c
c *** Background model grid data.
c
      real*4    prbg(nx_bg,ny_bg,nz_bg),     !Pressure (mb)
     .          htbg(nx_bg,ny_bg,nz_bg),     !Height (m)
     .          tpbg(nx_bg,ny_bg,nz_bg),     !Temperature (K)
     .          shbg(nx_bg,ny_bg,nz_bg),     !Specific humidity (kg/kg)
     .          uwbg(nx_bg,ny_bg,nz_bg),     !U-wind (m/s)
     .          vwbg(nx_bg,ny_bg,nz_bg),     !V-wind (m/s)
     .          mslpbg(nx_bg,ny_bg),         !mslp  (mb)
     .          wwbg(nx_bg,ny_bg,nz_bg),     !W-wind (pa/s)
     .          htbg_sfc(nx_bg,ny_bg),
     .          prbg_sfc(nx_bg,ny_bg),
     .          shbg_sfc(nx_bg,ny_bg),
     .          uwbg_sfc(nx_bg,ny_bg),
     .          vwbg_sfc(nx_bg,ny_bg),
     .          tpbg_sfc(nx_bg,ny_bg)

      real*4    lon0,lat1,lat2

      character*(*)   gproj,cmodel
      character*(*)   fullname
      character*(*)   bgpath
      character*(*)   fname_bg
      character*(*)   af_bg
      character*6     c6_maproj
      character*13    fname13
      character*13    fname9_to_wfo_fname13

      if (bgmodel .eq. 1) then     ! Process 60 km RUC data
          call read_ruc60_native(bgpath,fname_bg,af_bg,
     .               nx_bg,ny_bg,nz_bg,prbg,htbg,tpbg,shbg,uwbg,vwbg,
     .               gproj,lon0,istatus)

      elseif (bgmodel .eq. 2) then ! Process 48 km ETA conus-c grid data
          call read_eta_conusC(fullname,nx_bg,ny_bg,nz_bg,
     .                         htbg, prbg,tpbg,uwbg,vwbg,shbg,wwbg,
     .                         htbg_sfc,prbg_sfc,shbg_sfc,tpbg_sfc,
     .                         uwbg_sfc,vwbg_sfc,mslpbg,
     .                         istatus)

          if(istatus.eq.0) then
             call lprep_eta_conusc(nx_bg,ny_bg,nz_bg,prbg,tpbg,shbg
     +            ,tpbg_sfc,prbg_sfc,shbg_sfc
     +            ,gproj,lon0,lat1,lat2,istatus)
          endif
c
      elseif (bgmodel .eq. 4) then ! Process SBN Conus 211 data (Eta or RUC)

          ntbg=10
          fname13 = fname9_to_wfo_fname13(fname_bg(1:9))
          call get_sbn_dims(bgpath,fname13,nxbg,nybg,nzbg,ntbg)

          print*,'entering read_conus_211'
          call read_conus_211(bgpath,fname_bg(1:9),af_bg,
     .         nx_bg,ny_bg,nz_bg, nxbg,nybg,nzbg,ntbg,
     .         prbg,htbg,tpbg,shbg,uwbg,vwbg,wwbg,
     .         prbg_sfc,uwbg_sfc,vwbg_sfc,shbg_sfc,tpbg_sfc,
     .         mslpbg,gproj,lon0,lat1,lat2,1,istatus)
c
      elseif (bgmodel .eq. 5) then ! Process 40 km RUC data

          call read_ruc2_hybb(fullname,nx_bg,ny_bg,nz_bg,
     +         mslpbg,htbg,prbg,shbg,uwbg,vwbg,tpbg,wwbg,istatus)
          if(istatus.eq.0) then
             print*,'Read complete: entering prep'
             call lprep_ruc2_hybrid(nx_bg,ny_bg,nz_bg,htbg,prbg,shbg,
     +            uwbg,vwbg,tpbg,uwbg_sfc,vwbg_sfc,tpbg_sfc,prbg_sfc,
     +            shbg_sfc,htbg_sfc,gproj,lon0,lat1,lat2)
             print*,'Data prep complete'
          else
             print*,'lprep_ruc not called'
          endif
c
c ETA grib ingest currently disabled (J. Smart 9-4-98)
c Also, NOGAPS 2.5 degree obsolete.
c bgmodel 3 = FA (Taiwan). bgmodel 6 = NOGAPS1.0. bgmodel 8 = AVN 1.0 deg
c
      elseif (bgmodel .eq. 3 .or.
     .        bgmodel .eq. 6 .or.
     .        bgmodel .eq. 8) then ! Process AVN or NOGAPS1.0 grib data

             call read_dgprep(bgmodel,cmodel,bgpath
     .                 ,fname_bg,af_bg,nx_bg,ny_bg,nz_bg
     .                 ,prbg,htbg,tpbg,shbg,uwbg,vwbg,wwbg
     .                 ,htbg_sfc,prbg_sfc,shbg_sfc,tpbg_sfc
     .                 ,uwbg_sfc,vwbg_sfc,mslpbg
     .                 ,gproj,lon0,lat1,lat2,istatus)

      elseif (bgmodel .eq. 9) then ! Process NWS Conus data (RUC,ETA,NGM,AVN)
             call read_conus_nws(bgpath,fname_bg,af_bg,
     .               nx_bg,ny_bg,nz_bg,prbg,htbg,tpbg,shbg,uwbg,vwbg,
     .               gproj,lon0,lat1,lat2,istatus)
c
      endif

      if(istatus.ne. 0)then
         print*,'Error reading background model data in read_bgdata'
      endif

      return
      end
