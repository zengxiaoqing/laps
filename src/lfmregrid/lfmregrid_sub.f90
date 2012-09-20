
subroutine lfmregrid_sub(nx_bg,ny_bg,nzbg,fname_in,NX_L,NY_L,gproj,lfmprd_dir,laps_data_root,laps_reftime,laps_valtime)

use mem_namelist

character*256 fname_in, fullname_in, lfmprd_dir, laps_data_root
character*132 cmodel
character*2   gproj
character*256 syscmd

integer bgmodel

integer ct,n2df,n3df
parameter (n2df = 8)
parameter (n3df = 0)

real sgrid(NX_L,NY_L,n2df)
character*132 com2d(n2df)  
character*3 name2d(n2df)  

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
real usf_laps(NX_L,NY_L)
real vsf_laps(NX_L,NY_L)
real swi_laps(NX_L,NY_L)
real s8a_laps(NX_L,NY_L)
real rto_laps(NX_L,NY_L)

logical wrapped

fullname_in = trim(fname_in)//'.fsf'

call read_fuafsf_cdf(fname_in & 
                    ,nx_bg, ny_bg, nzbg &
                    ,htbg, pr, wwbg, shbg, tpbg, uwbg, vwbg &
                    ,uwbg_sfc, vwbg_sfc, tpbg_sfc, tdbg_sfc &     
                    ,prbg_sfc, mslpbg, htbg_sfc &
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

else

    write(6,*)' Calling bilinear_laps_2d for lmr'
    call bilinear_laps_2d(grx,gry,nx_bg,ny_bg,NX_L,NY_L,lmr,lmr_laps)

    write(6,*)' Calling bilinear_laps_2d for tsf'
    call bilinear_laps_2d(grx,gry,nx_bg,ny_bg,NX_L,NY_L,tpbg_sfc,tsf_laps)

    write(6,*)' Calling bilinear_laps_2d for dsf'
    call bilinear_laps_2d(grx,gry,nx_bg,ny_bg,NX_L,NY_L,tdbg_sfc,dsf_laps)

    write(6,*)' Calling bilinear_laps_2d for usf'
    call bilinear_laps_2d(grx,gry,nx_bg,ny_bg,NX_L,NY_L,uwbg_sfc,usf_laps)

    write(6,*)' Calling bilinear_laps_2d for vsf'
    call bilinear_laps_2d(grx,gry,nx_bg,ny_bg,NX_L,NY_L,vwbg_sfc,vsf_laps)

    write(6,*)' Calling bilinear_laps_2d for swi'
    call bilinear_laps_2d(grx,gry,nx_bg,ny_bg,NX_L,NY_L,swi,swi_laps)

    write(6,*)' Calling bilinear_laps_2d for rto'
    call bilinear_laps_2d(grx,gry,nx_bg,ny_bg,NX_L,NY_L,pcpbg,rto_laps)


endif

s8a_laps = r_missing_data
where (rto_laps(:,:) .ne. r_missing_data)
    rto_laps(:,:) = rto_laps(:,:) / 1000.
endwhere

write(6,*)' lmr ranges: ',minval(lmr),maxval(lmr),minval(lmr_laps),maxval(lmr_laps)
write(6,*)' tsf ranges: ',minval(tpbg_sfc),maxval(tpbg_sfc),minval(tsf_laps),maxval(tsf_laps)
write(6,*)' dsf ranges: ',minval(tdbg_sfc),maxval(tdbg_sfc),minval(dsf_laps),maxval(dsf_laps)
write(6,*)' usf ranges: ',minval(uwbg_sfc),maxval(uwbg_sfc),minval(usf_laps),maxval(usf_laps)
write(6,*)' vsf ranges: ',minval(vwbg_sfc),maxval(vwbg_sfc),minval(vsf_laps),maxval(vsf_laps)
write(6,*)' swi ranges: ',minval(swi),maxval(swi),minval(swi_laps),maxval(swi_laps)
write(6,*)' s8a ranges: ',minval(s8a),maxval(s8a),minval(s8a_laps),maxval(s8a_laps)
write(6,*)' rto ranges: ',minval(pcpbg),maxval(pcpbg),minval(rto_laps),maxval(rto_laps)

write(6,*)' Setup output arrays'
ct=1
sgrid(1:NX_L,1:NY_L,ct) = lmr_laps; name2d(ct)='LMR'; com2d(ct)='Composite Reflectivity'     ; ct=ct+1
sgrid(1:NX_L,1:NY_L,ct) = tsf_laps; name2d(ct)='TSF'; com2d(ct)='Sfc Temperature'            ; ct=ct+1
sgrid(1:NX_L,1:NY_L,ct) = dsf_laps; name2d(ct)='DSF'; com2d(ct)='Sfc Dewpoint Temperature'   ; ct=ct+1
sgrid(1:NX_L,1:NY_L,ct) = usf_laps; name2d(ct)='USF'; com2d(ct)='Sfc U-Component Wind'       ; ct=ct+1
sgrid(1:NX_L,1:NY_L,ct) = vsf_laps; name2d(ct)='VSF'; com2d(ct)='Sfc V-Component Wind'       ; ct=ct+1
sgrid(1:NX_L,1:NY_L,ct) = swi_laps; name2d(ct)='SWI'; com2d(ct)='Incoming SW Radiation'      ; ct=ct+1
sgrid(1:NX_L,1:NY_L,ct) = s8a_laps; name2d(ct)='S8A'; com2d(ct)='11u Brightness Temperature' ; ct=ct+1
sgrid(1:NX_L,1:NY_L,ct) = rto_laps; name2d(ct)='RTO'; com2d(ct)='Run Total Precip Accum'     ; ct=ct+1

if(ct-1 .ne. n2df)then
    write(6,*)' ERROR in lfmregrid_sub: value of ct is inconsistent with n2df ',ct,n2df
endif

NZ_L = 1
syscmd = 'rm -f '//trim(fullname_in)
call system (syscmd)                  

write(6,*)
write(6,*)' Calling output_laps_rg: NX_L,NY_L is ',NX_L,NY_L              
write(6,*)' Output directory is ',lfmprd_dir
call output_laps_rg(lfmprd_dir,laps_data_root,domnum_in,laps_reftime,laps_valtime,sgrid,name2d,com2d,n2df,n3df,NX_L,NY_L,NZ_L)

write(6,*)' End of lfmregrid_sub'

return

end

