
!GFORTRAN subroutine lfmregrid_sub(nx_bg,ny_bg,nzbg,nzbg_tp,nzbg_sh,nzbg_uv,nzbg_ww,fname_in,NX_L,NY_L,NZ_L,gproj,bgmodel,cmodel,laps_data_root,mtype,laps_reftime,laps_valtime,l_process_grib,l_process_cdf,l_grib_fua,l_grib_fsf,l_cdf_fua,l_cdf_fsf)
!GFORTRAN modifications begin
subroutine lfmregrid_sub(nx_bg,ny_bg,nzbg,nzbg_tp,nzbg_sh,nzbg_uv,nzbg_ww,&
                         fname_in,NX_L,NY_L,NZ_L,gproj,bgmodel,cmodel,&
                         laps_data_root,mtype,laps_reftime,laps_valtime,&
                         l_process_grib,l_process_cdf,l_grib_fua,l_grib_fsf,&
                         l_cdf_fua,l_cdf_fsf)
!GFORTRAN modifications end

!use mem_namelist
use storage_module, ONLY: get_plvls
integer, parameter :: maxbglvl = 52 ! Dimension is maxbglvl in 'degrib_nav' routine

character*256 fname_in, fullname_in, laps_data_root, fname_bg, bgpath
character*132 cmodel
character*4 af_bg
character*2   gproj
character*256 syscmd
character*30 mtype 

integer bgmodel

integer ct,n2df,n3df
parameter (n2df = 11)
parameter (n3df = 4)

real sgrid(NX_L,NY_L,n2df)
real pgrid(NX_L,NY_L,n3df*NZ_L)
character*132 com2d(n2df)  
character*132 com3d(n3df*NZ_L)  
character*3 name2d(n2df)  
character*3 name3d(n3df*NZ_L)  
integer lvls3d(n3df*NZ_L)  ! mb for pressure grid

real pres_1d(NZ_L) ! pa
real pr1d_pa(NZ_L) ! pa
real pr1d_mb(NZ_L) ! mb

! *** Background model grid data.
!
real  prbght(nx_bg,ny_bg,nzbg) !Pressure (mb) ht and temp
real  prbgsh(nx_bg,ny_bg,nzbg) !Pressure (mb) q
real  prbguv(nx_bg,ny_bg,nzbg) !Pressure (mb) u- v-components
real  prbgww(nx_bg,ny_bg,nzbg) !Pressure (mb) omega
real  pr(nzbg)

real  htbg(nx_bg,ny_bg,nzbg)   !Height (m)
real  tpbg(nx_bg,ny_bg,nzbg)   !Temperature (K)
real  shbg(nx_bg,ny_bg,nzbg)   !Specific humidity (kg/kg)
real  tdbg(nx_bg,ny_bg,nzbg)   !Dewpoint (K)
real  uwbg(nx_bg,ny_bg,nzbg)   !U-wind (m/s)
real  vwbg(nx_bg,ny_bg,nzbg)   !V-wind (m/s)
real  wwbg(nx_bg,ny_bg,nzbg)   !W-wind (pa/s)
real plvl_grib(maxbglvl) ! Dimension is maxbglvl in 'degrib_nav' routine

real  mslpbg(nx_bg,ny_bg)         !mslp  (mb)
real  htbg_sfc(nx_bg,ny_bg)
real  prbg_sfc(nx_bg,ny_bg)
real  shbg_sfc(nx_bg,ny_bg)       !Specific humidity (kg/kg)
real  uwbg_sfc(nx_bg,ny_bg)
real  vwbg_sfc(nx_bg,ny_bg)
real  tdbg_sfc(nx_bg,ny_bg)
real  tpbg_sfc(nx_bg,ny_bg)
real  t_at_sfc(nx_bg,ny_bg)
real  pcpbg(nx_bg,ny_bg)          !Precip at surface, ACPC (k/m^2)

!     Local input model variables for the time being
real r01(nx_bg,ny_bg)
real lmr(nx_bg,ny_bg)
real llr(nx_bg,ny_bg)
real swi(nx_bg,ny_bg)
real s8a(nx_bg,ny_bg)
real tpw(nx_bg,ny_bg)

real htvi(nx_bg,ny_bg,NZ_L)
real tpvi(nx_bg,ny_bg,NZ_L)
real shvi(nx_bg,ny_bg,NZ_L)
real uwvi(nx_bg,ny_bg,NZ_L)
real vwvi(nx_bg,ny_bg,NZ_L)
real wwvi(nx_bg,ny_bg,NZ_L)

real lat(NX_L,NY_L)         ! LAPS lat
real lon(NX_L,NY_L)         ! LAPS lon
real topo(NX_L,NY_L)        ! LAPS lon
real grx(NX_L,NY_L)         ! hinterp factor (background grid x index)
real gry(NX_L,NY_L)         ! hinterp factor (background grid y index)
real lmr_laps(NX_L,NY_L)
real tsf_laps(NX_L,NY_L)
real dsf_laps(NX_L,NY_L)
real psf_laps(NX_L,NY_L)
real usf_laps(NX_L,NY_L)
real vsf_laps(NX_L,NY_L)
real swi_laps(NX_L,NY_L)
real s8a_laps(NX_L,NY_L)
real tpw_laps(NX_L,NY_L)
real rto_laps(NX_L,NY_L)
real r01_laps(NX_L,NY_L)

real ssh, k_to_c

logical wrapped, l_process_grib, l_process_cdf
logical     l_grib_fua,l_grib_fsf,l_cdf_fua,l_cdf_fsf

write(6,*)                              
write(6,*)' Subroutine lfmregrid_sub...'
write(6,*)' nx_bg/ny_bg = ',nx_bg,ny_bg
write(6,*)' cmodel = ',cmodel          

write(6,*)' l_grib_fua,l_grib_fsf,l_cdf_fua,l_cdf_fsf:', l_grib_fua,l_grib_fsf,l_cdf_fua,l_cdf_fsf

call get_laps_domain(NX_L,NY_L,'nest7grid' &
                    ,lat,lon,topo,istatus)
if (istatus.lt.1)then
    print *,'Error reading lat, lon, topo data from get_laps_domain'
    stop
endif

call get_pres_1d(i4_valtime,NZ_L,pr1d_pa,istatus)
if(istatus.ne.1)then
   print*,'Error returned from get_pres_1d'
   print*,'Check pressures.nl or nk_laps in nest7grid.parms'
   stop
endif
pr1d_mb(:)=pr1d_pa(:)/100.  ! Pa to mb

if(l_process_grib .eqv. .true.)then
!   bgmodel = 13
!   if(mtype .eq. 'nam')then
!       cmodel = 'NAM'
!   else
!       cmodel = 'HRRR'
!   endif
!   write(6,*)' calling get_bkgd_mdl_info for cmodel: ',trim(cmodel)
!   call get_bkgd_mdl_info(bgmodel,cmodel,fname_in  &
!     ,nx_bg,ny_bg,nzbg,nzbg_tp,nzbg_sh,nzbg_uv,nzbg_ww &
!     ,gproj,dlat,dlon,centrallat,centrallon,dxbg,dybg &
!     ,Lat0,Lat1,Lon0,sw,ne,cgrddef,istatus)

    if(.true.)then
        call get_basename(fname_in,fname_bg)
        call get_directory_length(fname_in,lend)
        bgpath=fname_in(1:lend)
        af_bg=fname_bg(10:13)
        write(6,*)' calling read_bgdata: dims are ',nx_bg,ny_bg,nzbg,nzbg_tp,nzbg_sh,nzbg_uv,nzbg_ww          
        call read_bgdata(nx_bg,ny_bg, &
           nzbg,nzbg_tp,nzbg_sh,nzbg_uv,nzbg_ww,'lapsb' &
          ,bgpath,fname_bg,af_bg,fname_in,cmodel,bgmodel &
          ,prbght,prbgsh,prbguv,prbgww &
          ,htbg, tpbg,uwbg,vwbg,shbg,wwbg &
          ,htbg_sfc,prbg_sfc,shbg_sfc,tdbg_sfc,tpbg_sfc &
          ,t_at_sfc,uwbg_sfc,vwbg_sfc,mslpbg,pcpbg,lmr,tpw,swi,istatus)

        write(6,*)' returned from read_bgdata'

        istatus=ishow_timer()

        nx_pr = 1
        ny_pr = 1

        call get_plvls(plvl_grib, 100, nlvl_grib)
        write(6,*)' grib plvls info: ',nlvl_grib,plvl_grib(1:nlvl_grib)
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

        ixmin = 1
        ixmax = nx_bg
        iymin = 1
        iymax = ny_bg

        nz_laps = NZ_L
        write(6,*)' calling vinterp: dims are ',nx_bg,ny_bg,nzbg_ht,nz_laps              
        call vinterp(nz_laps,nx_bg,ny_bg,nx_pr,ny_pr &
               ,ixmin,ixmax,iymin,iymax &
               ,nzbg_ht,nzbg_tp,nzbg_sh,nzbg_uv,nzbg_ww &
               ,pr1d_mb,prbght,prbgsh,prbguv,prbgww &
               ,htbg,tpbg,shbg,uwbg,vwbg,wwbg &
               ,htvi,tpvi,shvi,uwvi,vwvi,wwvi)

        nx_laps = NX_L
        ny_laps = NY_L

!       if((l_cdf_fua .eqv. .true.) .OR. (l_grib_fua .eqv. .true.))then
        if(.false.)then

        write(6,*)' calling init_hinterp'
        wrapped = .false.
        bgmodel = 0
        call init_hinterp(nx_bg,ny_bg,NX_L,NY_L,gproj, &
                          lat,lon,grx,gry,bgmodel,cmodel,wrapped)

        print*,'use bilinear_laps_3d for T ',trim(cmodel)
        call bilinear_laps_3df(grx,gry,nx_bg,ny_bg &
                              ,nx_laps,ny_laps,nz_laps,tpvi,tp)

        istatus=ishow_timer()
        print*,'use bilinear_laps_3d for Q ',trim(cmodel)
        call bilinear_laps_3df(grx,gry,nx_bg,ny_bg &
                              ,nx_laps,ny_laps,nz_laps,shvi,sh)

        endif

        if(istatus .eq. 1)then
            write(6,*)' calling write_lga'            
            call write_lga(nx_laps,ny_laps,nz_laps,laps_reftime, &
                           bgvalid,cmodel,missingflag,pr1d_mb,ht,tp,sh,uw,vw,ww,istatus)
            if(istatus.ne.1)then
             print*,'Error writing lga - returning to main'
            endif

        endif
    endif
endif

if(l_process_cdf .eqv. .true.)then
!   fullname_in = trim(fname_in)//'.fsf'

    write(6,*)' Calling read_fuafsf_cdf'

!   Initialize variables
    tpw = r_missing_data

    call read_fuafsf_cdf(fname_in & 
                    ,nx_bg, ny_bg, nzbg &
                    ,htbg, pr, wwbg, shbg, tpbg, uwbg, vwbg &
                    ,uwbg_sfc, vwbg_sfc, tpbg_sfc, tdbg_sfc &     
                    ,prbg_sfc, mslpbg, htbg_sfc &
                    ,r01 &
                    ,pcpbg &
                    ,lmr, llr, s8a, swi, tpw &
                    ,istatus) 
    if(istatus.ne.1 .AND. istatus.ne.-1)then
        print*,'Error returned: read_fuafsf_cdf'
        return
    endif

endif

! if(l_process_cdf .eqv. .true.)then
if(.true.)then
    wrapped = .false.
    bgmodel = 0

    write(6,*)' Calling init_hinterp'

    call init_hinterp(nx_bg,ny_bg,NX_L,NY_L,gproj, &
                  lat,lon,grx,gry,bgmodel,cmodel,wrapped)

    if((l_grib_fua .eqv. .true.) .OR. (l_cdf_fua .eqv. .true.))then
        call get_pres_1d(i4_valtime,NZ_L,pres_1d,istatus)
        lz = NZ_L
        ct=1

        print*,'use bilinear_laps_3df for T3 starting at pgrid level',ct                        
        call bilinear_laps_3df(grx,gry,nx_bg,ny_bg &
                              ,NX_L,NY_L,NZ_L,tpbg,pgrid(1,1,ct))
        name3d(ct:ct+lz-1)='T3 '; com3d(ct:ct+lz-1)='Temperature'; lvls3d(ct:ct+lz-1)=nint(pres_1d(lz:1:-1)/100.); ct=ct+lz   

        if(trim(cmodel) .eq. 'HRRR')then
!         Convert sh from dpt into actual sh
          itesth = 1; jtesth = 1
          itestl = 506; jtestl = 128
          if(NX_L .gt. itestl .and. NY_L .gt. jtestl)then ! DFW LAPS grid location (HWT domain)
              itesth = nint(grx(506,128))
              jtesth = nint(gry(506,128))
              write(6,*)' Convert sh from dpt to actual sh at lat/lon ',lat(itestl,jtestl),lon(itestl,jtestl)            
              write(6,*)' LAPS gridpt',itestl,jtestl,' HRRR gridpt',itesth,jtesth
          endif
              
          do k = 1,NZ_L
            kflip = (NZ_L+1) - k
            p_mb = pres_1d(kflip) / 100.
            do i = 1,nx_bg
            do j = 1,ny_bg
                sh_orig = shbg(i,j,k)
                td_c = k_to_c(shbg(i,j,k))
                shbg(i,j,k) = ssh(p_mb,td_c) / 1000.
                if(i .eq. itesth .and. j .eq. jtesth)then
                    write(6,11)p_mb,tpbg(i,j,k),sh_orig,td_c,shbg(i,j,k)
                endif
11              format(' p,t,sh_orig,td_c,sh: ',4f10.3,f10.4)
            enddo ! j 
            enddo ! i
          enddo ! k
        endif

        istatus=ishow_timer()

        print*,'use bilinear_laps_3df for SH starting at pgrid level',ct                        
        write(6,*)shbg(1,1,:)
        call bilinear_laps_3df(grx,gry,nx_bg,ny_bg &
                              ,NX_L,NY_L,NZ_L,shbg,pgrid(1,1,ct))
        write(6,*)pgrid(1,1,ct:ct+lz-1)
        name3d(ct:ct+lz-1)='SH '; com3d(ct:ct+lz-1)='Specific Humidity'; lvls3d(ct:ct+lz-1)=nint(pres_1d(lz:1:-1)/100.); ct=ct+lz   

        istatus=ishow_timer()

        print*,'use bilinear_laps_3df for U3 starting at pgrid level',ct                        
        call bilinear_laps_3df(grx,gry,nx_bg,ny_bg &
                              ,NX_L,NY_L,NZ_L,uwbg,pgrid(1,1,ct))
        name3d(ct:ct+lz-1)='U3 '; com3d(ct:ct+lz-1)='U-component Wind'; lvls3d(ct:ct+lz-1)=nint(pres_1d(lz:1:-1)/100.); ct=ct+lz   

        istatus=ishow_timer()

        print*,'use bilinear_laps_3df for V3 starting at pgrid level',ct                        
        call bilinear_laps_3df(grx,gry,nx_bg,ny_bg &
                              ,NX_L,NY_L,NZ_L,vwbg,pgrid(1,1,ct))
        name3d(ct:ct+lz-1)='V3 '; com3d(ct:ct+lz-1)='V-component Wind'; lvls3d(ct:ct+lz-1)=nint(pres_1d(lz:1:-1)/100.); ct=ct+lz   

        istatus=ishow_timer()

    else
        write(6,*)' Skip 3D interpolation - FUA processing set to FALSE'

    endif

    if(.false.)then
  
        write(6,*)' Calling hinterp_field_2d for lmr'
        call hinterp_field_2d(nx_bg,ny_bg,NX_L,NY_L,1, &
                              grx,gry,lmr,lmr_laps,wrapped)

        write(6,*)' Calling hinterp_field_2d for tsf'
        call hinterp_field_2d(nx_bg,ny_bg,NX_L,NY_L,1, &
                              grx,gry,tpbg_sfc,tsf_laps,wrapped)

        write(6,*)' Calling hinterp_field_2d for dsf'
        call hinterp_field_2d(nx_bg,ny_bg,NX_L,NY_L,1, &
                              grx,gry,tdbg_sfc,dsf_laps,wrapped)

        write(6,*)' Calling hinterp_field_2d for psf'
        call hinterp_field_2d(nx_bg,ny_bg,NX_L,NY_L,1, &
                              grx,gry,prbg_sfc,psf_laps,wrapped)

        write(6,*)' Calling hinterp_field_2d for usf'
        call hinterp_field_2d(nx_bg,ny_bg,NX_L,NY_L,1, &
                              grx,gry,uwbg_sfc,usf_laps,wrapped)

        write(6,*)' Calling hinterp_field_2d for vsf'
        call hinterp_field_2d(nx_bg,ny_bg,NX_L,NY_L,1, &
                              grx,gry,vwbg_sfc,vsf_laps,wrapped)

        write(6,*)' Calling hinterp_field_2d for swi'
        call hinterp_field_2d(nx_bg,ny_bg,NX_L,NY_L,1, &
                             grx,gry,swi,swi_laps,wrapped)

        write(6,*)' Calling hinterp_field_2d for rto'
        call hinterp_field_2d(nx_bg,ny_bg,NX_L,NY_L,1, &
                             grx,gry,pcpbg,rto_laps,wrapped)

        write(6,*)' Calling hinterp_field_2d for r01'
        call hinterp_field_2d(nx_bg,ny_bg,NX_L,NY_L,1, &
                             grx,gry,pcpbg,r01_laps,wrapped)

        write(6,*)' Calling hinterp_field_2d for tpw'
        call hinterp_field_2d(nx_bg,ny_bg,NX_L,NY_L,1, &
                             grx,gry,tpw,tpw_laps,wrapped)

    else

        write(6,*)' Calling bilinear_laps_2d for lmr'
        call bilinear_laps_2d(grx,gry,nx_bg,ny_bg,NX_L,NY_L,lmr,lmr_laps)

        write(6,*)' Calling bilinear_laps_2d for tsf'
        call bilinear_laps_2d(grx,gry,nx_bg,ny_bg,NX_L,NY_L,tpbg_sfc,tsf_laps)

        write(6,*)' Calling bilinear_laps_2d for dsf'
        call bilinear_laps_2d(grx,gry,nx_bg,ny_bg,NX_L,NY_L,tdbg_sfc,dsf_laps)

        write(6,*)' Calling bilinear_laps_2d for psf'
        call bilinear_laps_2d(grx,gry,nx_bg,ny_bg,NX_L,NY_L,prbg_sfc,psf_laps)

        write(6,*)' Calling bilinear_laps_2d for usf'
        call bilinear_laps_2d(grx,gry,nx_bg,ny_bg,NX_L,NY_L,uwbg_sfc,usf_laps)

        write(6,*)' Calling bilinear_laps_2d for vsf'
        call bilinear_laps_2d(grx,gry,nx_bg,ny_bg,NX_L,NY_L,vwbg_sfc,vsf_laps)

        write(6,*)' Calling bilinear_laps_2d for swi'
        call bilinear_laps_2d(grx,gry,nx_bg,ny_bg,NX_L,NY_L,swi,swi_laps)

        write(6,*)' Calling bilinear_laps_2d for rto'
        call bilinear_laps_2d(grx,gry,nx_bg,ny_bg,NX_L,NY_L,pcpbg,rto_laps)

        write(6,*)' Calling bilinear_laps_2d for r01'
        call bilinear_laps_2d(grx,gry,nx_bg,ny_bg,NX_L,NY_L,pcpbg,r01_laps)

        write(6,*)' Calling bilinear_laps_2d for tpw'
        call bilinear_laps_2d(grx,gry,nx_bg,ny_bg,NX_L,NY_L,tpw,tpw_laps)


    endif

    s8a_laps = r_missing_data

    where (rto_laps(:,:) .ne. r_missing_data)
        rto_laps(:,:) = rto_laps(:,:) / 1000.
    endwhere
    where (r01_laps(:,:) .ne. r_missing_data)
        r01_laps(:,:) = r01_laps(:,:) / 1000.
    endwhere
    where (tpw_laps(:,:) .ne. r_missing_data)
        tpw_laps(:,:) = tpw_laps(:,:) / 1000.
    endwhere

    write(6,*)' lmr ranges: ',minval(lmr),maxval(lmr),minval(lmr_laps),maxval(lmr_laps)
    write(6,*)' tsf ranges: ',minval(tpbg_sfc),maxval(tpbg_sfc),minval(tsf_laps),maxval(tsf_laps)
    write(6,*)' dsf ranges: ',minval(tdbg_sfc),maxval(tdbg_sfc),minval(dsf_laps),maxval(dsf_laps)
    write(6,*)' psf ranges: ',minval(prbg_sfc),maxval(prbg_sfc),minval(psf_laps),maxval(psf_laps)
    write(6,*)' usf ranges: ',minval(uwbg_sfc),maxval(uwbg_sfc),minval(usf_laps),maxval(usf_laps)
    write(6,*)' vsf ranges: ',minval(vwbg_sfc),maxval(vwbg_sfc),minval(vsf_laps),maxval(vsf_laps)
    write(6,*)' swi ranges: ',minval(swi),maxval(swi),minval(swi_laps),maxval(swi_laps)
    write(6,*)' s8a ranges: ',minval(s8a),maxval(s8a),minval(s8a_laps),maxval(s8a_laps)
    write(6,*)' rto ranges: ',minval(pcpbg),maxval(pcpbg),minval(rto_laps),maxval(rto_laps)
    write(6,*)' r01 ranges: ',minval(r01),maxval(r01),minval(r01_laps),maxval(r01_laps)
    write(6,*)' tpw ranges: ',minval(tpw),maxval(tpw),minval(tpw_laps),maxval(tpw_laps)

    write(6,*)' Setup output arrays'
    ct=1
    sgrid(1:NX_L,1:NY_L,ct) = lmr_laps; name2d(ct)='LMR'; com2d(ct)='Composite Reflectivity'     ; ct=ct+1
    sgrid(1:NX_L,1:NY_L,ct) = tsf_laps; name2d(ct)='TSF'; com2d(ct)='Sfc Temperature'            ; ct=ct+1
    sgrid(1:NX_L,1:NY_L,ct) = dsf_laps; name2d(ct)='DSF'; com2d(ct)='Sfc Dewpoint Temperature'   ; ct=ct+1
    sgrid(1:NX_L,1:NY_L,ct) = psf_laps; name2d(ct)='PSF'; com2d(ct)='Surface Pressure'           ; ct=ct+1
    sgrid(1:NX_L,1:NY_L,ct) = usf_laps; name2d(ct)='USF'; com2d(ct)='Sfc U-Component Wind'       ; ct=ct+1
    sgrid(1:NX_L,1:NY_L,ct) = vsf_laps; name2d(ct)='VSF'; com2d(ct)='Sfc V-Component Wind'       ; ct=ct+1
    sgrid(1:NX_L,1:NY_L,ct) = swi_laps; name2d(ct)='SWI'; com2d(ct)='Incoming SW Radiation'      ; ct=ct+1
    sgrid(1:NX_L,1:NY_L,ct) = s8a_laps; name2d(ct)='S8A'; com2d(ct)='11u Brightness Temperature' ; ct=ct+1
    sgrid(1:NX_L,1:NY_L,ct) = rto_laps; name2d(ct)='RTO'; com2d(ct)='Run Total Precip Accum'     ; ct=ct+1
    sgrid(1:NX_L,1:NY_L,ct) = r01_laps; name2d(ct)='R01'; com2d(ct)='Incremental Precip Accum'   ; ct=ct+1
    sgrid(1:NX_L,1:NY_L,ct) = tpw_laps; name2d(ct)='TPW'; com2d(ct)='Precipitable Water Depth'   ; ct=ct+1

    if(ct-1 .ne. n2df)then
        write(6,*)' ERROR in lfmregrid_sub: value of ct is inconsistent with n2df ',ct,n2df
    endif

!   if(n3df .eq. 0)NZ_L = 1
!   syscmd = 'rm -f '//trim(fullname_in)
!   call system (syscmd)                  

    write(6,*)
    write(6,*)' Calling output_laps_rg: NX_L,NY_L,NZ_L is ',NX_L,NY_L,NZ_L              
    write(6,*)' laps_data_root is: ',laps_data_root
    write(6,*)' pgrid is ',pgrid(1,1,:) 
!GFORTRAN    call output_laps_rg(laps_data_root,mtype,domnum_in,laps_reftime,laps_valtime,pgrid,sgrid,name2d,name3d,com2d,com3d,lvls3d,n2df,n3df,pres_1d,NX_L,NY_L,NZ_L)
!GFORTRAN modifications begin
    call output_laps_rg(laps_data_root,mtype,domnum_in,laps_reftime,&
                        laps_valtime,pgrid,sgrid,name2d,name3d,com2d,com3d,&
                        lvls3d,n2df,n3df,pres_1d,NX_L,NY_L,NZ_L)
!GFORTRAN modifications end

endif ! .true. (formerly l_process_cdf)

write(6,*)' End of lfmregrid_sub'

return

end

subroutine get_basename(fullname,basename)

character*(*) fullname,basename

call get_directory_length(fullname,lend)
call s_len(fullname,lenfn)

basename = fullname(lend+1:lenfn)

write(6,*)' basename = ',trim(basename)

return
end
