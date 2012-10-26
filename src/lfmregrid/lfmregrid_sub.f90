
subroutine lfmregrid_sub(nx_bg,ny_bg,nzbg,fname_in,NX_L,NY_L,NZ_L,gproj,laps_data_root,mtype,laps_reftime,laps_valtime,l_process_grib,l_process_cdf,l_grib_fua,l_grib_fsf,l_cdf_fua,l_cdf_fsf)

use mem_namelist

character*256 fname_in, fullname_in, laps_data_root, fname_bg, bgpath
character*132 cmodel
character*4 af_bg
character*2   gproj
character*256 syscmd
character*30 mtype 

integer bgmodel

integer ct,n2df,n3df
parameter (n2df = 11)
parameter (n3df = 2)

real sgrid(NX_L,NY_L,n2df)
real pgrid(NX_L,NY_L,n3df*NZ_L)
character*132 com2d(n2df)  
character*132 com3d(n3df*NZ_L)  
character*3 name2d(n2df)  
character*3 name3d(n3df*NZ_L)  
integer lvls3d(n3df*NZ_L)  ! mb for pressure grid

real pres_1d(NZ_L) ! pa

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

real lat(NX_L,NY_L)         ! LAPS lat
real lon(NX_L,NY_L)         ! LAPS lon
real topo(NX_L,NY_L)        ! LAPS lon
real grx(NX_L,NY_L)         ! hinterp factor
real gry(NX_L,NY_L)         ! hinterp factor
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

real make_ssh, k_to_c

logical wrapped, l_process_grib, l_process_cdf
logical     l_grib_fua,l_grib_fsf,l_cdf_fua,l_cdf_fsf

write(6,*)' Subroutine lfmregrid_sub'

if(l_process_grib .eqv. .true.)then
    bgmodel = 13
    cmodel = 'HRRR'
    write(6,*)' calling get_bkgd_mdl_info'
    call get_bkgd_mdl_info(bgmodel,cmodel,fname_in  &
      ,nx_bg,ny_bg,nzbg,nzbg_tp,nzbg_sh,nzbg_uv,nzbg_ww &
      ,gproj,dlat,dlon,centrallat,centrallon,dxbg,dybg &
      ,Lat0,Lat1,Lon0,sw,ne,cgrddef,istatus)

    if(istatus .eq. 1)then
        call get_basename(fname_in,fname_bg)
        call get_directory_length(fname_in,lend)
        bgpath=fname_in(1:lend)
        af_bg=fname_bg(10:13)
        write(6,*)' calling read_bgdata'          
        call read_bgdata(nx_bg,ny_bg, &
           nzbg,nzbg_tp,nzbg_sh,nzbg_uv,nzbg_ww,'lapsb' &
          ,bgpath,fname_bg,af_bg,fname_in,cmodel,bgmodel &
          ,prbght,prbgsh,prbguv,prbgww &
          ,htbg, tpbg,uwbg,vwbg,shbg,wwbg &
          ,htbg_sfc,prbg_sfc,shbg_sfc,tdbg_sfc,tpbg_sfc &
          ,t_at_sfc,uwbg_sfc,vwbg_sfc,mslpbg,pcpbg,istatus)

        istatus=ishow_timer()

        nx_pr = 1
        ny_pr = 1

        write(6,*)' calling vinterp'              
        call vinterp(nz_laps,nx_bg,ny_bg,nx_pr,ny_pr &
               ,ixmin,ixmax,iymin,iymax &
               ,nzbg_ht,nzbg_tp,nzbg_sh,nzbg_uv,nzbg_ww &
               ,pr1d_mb,prbght,prbgsh,prbguv,prbgww &
               ,htbg,tpbg,shbg,uwbg,vwbg,wwbg &
               ,htvi,tpvi,shvi,uwvi,vwvi,wwvi)

        print*,'use bilinear_laps_3d for T ',cmodel(1:ic)
        call bilinear_laps_3df(grx,gry,nx_bg,ny_bg &
                              ,nx_laps,ny_laps,nz_laps,tpvi,tp)

        istatus=ishow_timer()
        print*,'use bilinear_laps_3d for Q ',cmodel(1:ic)
        call bilinear_laps_3df(grx,gry,nx_bg,ny_bg &
                              ,nx_laps,ny_laps,nz_laps,shvi,sh)

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

    call get_laps_domain(NX_L,NY_L,'nest7grid' &
                        ,lat,lon,topo,istatus)
    if (istatus.lt.1)then
        print *,'Error reading lat, lon, topo data from get_laps_domain'
        stop
    endif

    wrapped = .false.
    bgmodel = 0

    write(6,*)' Calling init_hinterp'

    call init_hinterp(nx_bg,ny_bg,NX_L,NY_L,gproj, &
                  lat,lon,grx,gry,bgmodel,cmodel,wrapped)

    if(l_cdf_fua .eqv. .true.)then
        call get_pres_1d(i4_valtime,NZ_L,pres_1d,istatus)
        lz = NZ_L
        ct=1

        print*,'use bilinear_laps_3df for T3 starting at pgrid level',ct                        
        call bilinear_laps_3df(grx,gry,nx_bg,ny_bg &
                              ,NX_L,NY_L,NZ_L,tpbg,pgrid(1,1,ct))
        name3d(ct:ct+lz-1)='T3 '; com3d(ct:ct+lz-1)='Temperature'; lvls3d(ct:ct+lz-1)=nint(pres_1d(lz:1:-1)/100.); ct=ct+lz   

!       Convert sh from dpt into actual sh
        do k = 1,NZ_L
            p_mb = pres_1d(k) / 100.
            do i = 1,nx_bg
            do j = 1,ny_bg
                sh_orig = shbg(i,j,k)
                tamb_c = k_to_c(shbg(i,j,k))
                shbg(i,j,k) = make_ssh(p_mb,tamb_c,1.0,-132.) / 1000.
                if(i .eq. 1 .and. j .eq. 1)write(6,*)' sh_orig,tamb_c,sh: ',sh_orig,tamb_c,shbg(i,j,k)
            enddo ! j 
            enddo ! i
        enddo ! k

        print*,'use bilinear_laps_3df for SH starting at pgrid level',ct                        
        write(6,*)shbg(1,1,:)
        call bilinear_laps_3df(grx,gry,nx_bg,ny_bg &
                              ,NX_L,NY_L,NZ_L,shbg,pgrid(1,1,ct))
        write(6,*)pgrid(1,1,ct:ct+lz-1)
        name3d(ct:ct+lz-1)='SH '; com3d(ct:ct+lz-1)='Specific Humidity'; lvls3d(ct:ct+lz-1)=nint(pres_1d(lz:1:-1)/100.); ct=ct+lz   

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
    call output_laps_rg(laps_data_root,mtype,domnum_in,laps_reftime,laps_valtime,pgrid,sgrid,name2d,name3d,com2d,com3d,lvls3d,n2df,n3df,pres_1d,NX_L,NY_L,NZ_L)

endif ! l_process_cdf

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
