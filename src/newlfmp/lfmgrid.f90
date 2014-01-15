module lfmgrid

use map_utils

implicit none

save

integer, parameter :: maxdomains=10
real, parameter :: rmsg=1.e37

integer :: nx,ny,nz,lx,ly,lz,nvar2d,nvar3d,nvar3dh,nvar2dout,nvar3dout

type(proj_info) :: proj

logical :: out_cdf=.true.,out_grib=.false.,out_v5d=.true.
logical :: make_micro=.true.,make_firewx=.true.,make_points(maxdomains)=.false.
logical :: verbose=.true.
logical :: realtime=.true.
logical :: write_to_lapsdir=.true.
logical :: make_donefile=.true.
logical :: large_pgrid=.false.
logical :: large_ngrid=.false.
logical :: l_process_p=.false.
logical :: l_process_t=.false.
logical :: l_process_mr=.false.
logical :: l_process_z=.false.
logical :: l_process_uv=.false.
logical :: l_process_w=.false.

integer :: domnum,fcsttime,laps_reftime,laps_valtime,precip_dt=3600  ! in seconds
integer :: k_micro=999
integer :: n3d_pts_thr=33000000
integer :: p3d_pts_thr=22000000
integer :: advection_time = 0
integer :: i4_adv_cld = -1
integer :: i4_adv_pcp = -1

character(len=256) :: filename,filename0
character(len=32)  :: mtype
character(len=20)  :: c_m2z='rams'
character(len=3)   :: domnum_fstr
character(len=4)   :: wrf_version='3'


! Point forecast variables.

integer :: point_tz_utcoffset=0
character(len=5) :: point_vent_units='M^2/S'
character(len=3) :: point_tz_label='UTC'
character(len=3) :: point_windspd_units='KTS'
character(len=1) :: point_temp_units='F'

! Native grid map specifications.

real :: ngrid_spacingx,ngrid_spacingy,nstdlon,ntruelat1,ntruelat2
character(len=32) :: nprojection

! Laps grid map specifications.

real :: grid_spacing,stdlon,truelat1,truelat2
character(len=32) :: projection

!real, target, dimension(:,:,:) :: ngrid   ! Native model grid
!real, target, dimension(:,:,:) :: hgrid   ! Horizontally interpolated 3d grid
!real, target, dimension(:,:,:) :: pgrid   ! Isobaric output grid
!real, target, dimension(:,:,:) :: sgrid   ! Horizontally interpolated sfc output grid
real, target, allocatable, dimension(:,:,:) :: ngrid   ! Native model grid
real, target, allocatable, dimension(:,:,:) :: hgrid   ! Horizontally interpolated 3d grid
real, target, allocatable, dimension(:,:,:) :: pgrid   ! Isobaric output grid
real, target, allocatable, dimension(:,:,:) :: sgrid   ! Horizontally interpolated sfc output grid

! Native grid variables.

real, allocatable, dimension(:,:) :: nlat,nlon
real, pointer, dimension(:,:) ::  &
       nzsfc    ,npsfc    ,ntsfc    ,nmrsfc   ,nusfc     &
      ,nvsfc    ,nwsfc    ,nground_t,npblhgt  ,npcp_tot  &
      ,nlwout   ,nswout   ,nlwdown  ,nswdown  ,nshflux   &
      ,nlhflux  
real, pointer, dimension(:,:,:) :: &
       npsig    ,nzsig    ,ntsig    ,nmrsig   ,nusig     &
      ,nvsig    ,nwsig    ,ntkesig                       &
      ,ncldliqmr_sig,ncldicemr_sig,nrainmr_sig,nsnowmr_sig,ngraupelmr_sig &
      ,nrefl_sig

! Horizontally interpolated grid variables.

real, pointer, dimension(:,:,:) :: &
       hpsig    ,hzsig    ,htsig    ,hmrsig   ,husig     &
      ,hvsig    ,hwsig    ,htkesig                       &
      ,hcldliqmr_sig,hcldicemr_sig,hrainmr_sig,hsnowmr_sig,hgraupelmr_sig  &
      ,hrefl_sig,hzdr_sig,hldr_sig 

! Laps isobaric grid variables.

real, allocatable, dimension(:) :: lprs,lprsl
real, pointer, dimension(:,:,:) ::  &
       zprs         ,rhprs        ,tprs         ,shprs        ,uprs          &
      ,vprs         ,wprs         ,omprs        ,cldliqmr_prs  &
      ,cldicemr_prs ,rainmr_prs   ,snowmr_prs   ,graupelmr_prs,refl_prs      &
      ,zdr_prs      ,ldr_prs                                                 &
      ,pcptype_prs  

integer,            allocatable, dimension(:) :: lvls3d
character(len=3),   allocatable, dimension(:) :: name3d
character(len=10),  allocatable, dimension(:) :: units3d
character(len=4),   allocatable, dimension(:) :: lvltype3d
character(len=132), allocatable, dimension(:) :: com3d

! Laps surface grid variables.

real, allocatable, dimension(:,:) :: llat,llon
real, pointer, dimension(:,:) ::  &
       zsfc         ,psfc         ,tsfc         ,mrsfc        ,usfc        &
      ,vsfc         ,wsfc         ,ground_t     ,pblhgt       ,rhsfc       &  
      ,pcp_tot      ,thetasfc     ,ztw0         ,ztw1                      &
      ,thetaesfc    ,tdsfc        ,redp         ,pmsl         ,upbl        &
      ,vpbl         ,u80          ,v80          ,clwmrsfc     ,icemrsfc     ,snowmrsfc    ,rainmrsfc   &
      ,graupmrsfc   ,cldbase      ,cldtop       ,cldamt       ,ceiling     &
      ,intliqwater  ,intcldice    ,totpcpwater  ,max_refl     ,echo_tops    ,refl_sfc    &
      ,pcptype_sfc  ,pcp_inc      ,snow_inc     ,snow_tot &
      ,srhel        ,uhel         ,cape         ,cin          ,liftedind   &
      ,visibility   ,heatind      ,lwout        ,swout        ,lwdown      &
      ,swdown       ,shflux       ,lhflux       ,vnt_index    ,ham_index   &
      ,hah_index    ,fwi_index    ,fwx_index    ,upflux       ,bt11u       &
      ,cldalb       ,simvis
integer,            allocatable, dimension(:) :: lvls2d
character(len=3),   allocatable, dimension(:) :: name2d
character(len=10),  allocatable, dimension(:) :: units2d
character(len=4),   allocatable, dimension(:) :: lvltype2d
character(len=132), allocatable, dimension(:) :: com2d

! Grib parameters.

integer :: table_version=2,center_id=59,subcenter_id=2,process_id=255  &
          ,funit,startbyte,nbytes,igds(18)
integer, allocatable, dimension(:) :: param,leveltype,level1,level2  &
                                     ,timerange,scalep10             &
                                     ,paramua,scalep10ua
character(len=256) :: gribfile
logical, allocatable, dimension(:) :: gribit,gribitua

contains

!===============================================================================

subroutine alloc_native_grid

implicit none

integer :: ct

if (trim(mtype) /= 'st4') then
  nvar2d=16
  if(.not. large_ngrid)then
      nvar3d=7 ! add 1 for tkesig if needed 
  else ! large_ngrid
    if(.not. large_pgrid)then ! process U,V
      nvar3d=6
    else
      nvar3d=4
    endif
  endif
  if (make_micro) then
      nvar3d=nvar3d+5
      if(c_m2z == 'wrf') then
          nvar3d=nvar3d+1
      endif
  endif

  allocate(nlat(nx,ny),nlon(nx,ny))

  print*,'Native grid allocation 2d/3d/2dpts/3dpts = ',nvar2d,nvar3d,nx*ny*nvar2d,nx*ny*nz*nvar3d 

  allocate(ngrid(nx,ny,nvar2d+nvar3d*nz))

  ngrid=rmsg

  ct=1

  nzsfc    =>ngrid(1:nx,1:ny,ct); ct=ct+1
  npsfc    =>ngrid(1:nx,1:ny,ct); ct=ct+1
  ntsfc    =>ngrid(1:nx,1:ny,ct); ct=ct+1
  nmrsfc   =>ngrid(1:nx,1:ny,ct); ct=ct+1
  nusfc    =>ngrid(1:nx,1:ny,ct); ct=ct+1
  nvsfc    =>ngrid(1:nx,1:ny,ct); ct=ct+1
  nwsfc    =>ngrid(1:nx,1:ny,ct); ct=ct+1
  nground_t=>ngrid(1:nx,1:ny,ct); ct=ct+1
  npblhgt  =>ngrid(1:nx,1:ny,ct); ct=ct+1
  npcp_tot =>ngrid(1:nx,1:ny,ct); ct=ct+1
  nlwout   =>ngrid(1:nx,1:ny,ct); ct=ct+1
  nswout   =>ngrid(1:nx,1:ny,ct); ct=ct+1
  nlwdown  =>ngrid(1:nx,1:ny,ct); ct=ct+1
  nswdown  =>ngrid(1:nx,1:ny,ct); ct=ct+1
  nshflux  =>ngrid(1:nx,1:ny,ct); ct=ct+1
  nlhflux  =>ngrid(1:nx,1:ny,ct); ct=ct+1

  npsig    =>ngrid(1:nx,1:ny,ct:ct+nz-1); ct=ct+nz; l_process_p = .true.
  nzsig    =>ngrid(1:nx,1:ny,ct:ct+nz-1); ct=ct+nz; l_process_z = .true.
  ntsig    =>ngrid(1:nx,1:ny,ct:ct+nz-1); ct=ct+nz; l_process_t = .true.
  nmrsig   =>ngrid(1:nx,1:ny,ct:ct+nz-1); ct=ct+nz; l_process_mr = .true.

  if(.not. large_pgrid)then
      nusig    =>ngrid(1:nx,1:ny,ct:ct+nz-1); ct=ct+nz; l_process_uv = .true.
      nvsig    =>ngrid(1:nx,1:ny,ct:ct+nz-1); ct=ct+nz
      if(.not. large_ngrid)then   
          nwsig    =>ngrid(1:nx,1:ny,ct:ct+nz-1); ct=ct+nz; l_process_w = .true.
      endif
  endif

  if (make_micro) then
    k_micro = ct
    ncldliqmr_sig =>ngrid(1:nx,1:ny,ct:ct+nz-1); ct=ct+nz
    ncldicemr_sig =>ngrid(1:nx,1:ny,ct:ct+nz-1); ct=ct+nz
    nrainmr_sig   =>ngrid(1:nx,1:ny,ct:ct+nz-1); ct=ct+nz
    nsnowmr_sig   =>ngrid(1:nx,1:ny,ct:ct+nz-1); ct=ct+nz
    ngraupelmr_sig=>ngrid(1:nx,1:ny,ct:ct+nz-1); ct=ct+nz
    if(c_m2z == 'wrf') then
       nrefl_sig  =>ngrid(1:nx,1:ny,ct:ct+nz-1); ct=ct+nz
    endif
  endif
else
  nvar2d=1
  nvar3d=0

  allocate(nlat(nx,ny),nlon(nx,ny))

  allocate(ngrid(nx,ny,nvar2d+nvar3d*nz))

  ngrid=rmsg

  ct=1

  npcp_tot =>ngrid(1:nx,1:ny,ct); ct=ct+1
endif

return
end subroutine

!===============================================================================

subroutine alloc_hinterp_grid

implicit none

integer :: ct

if(.not. large_ngrid)then ! state variables
   nvar3dh=7
else
   if(.not. large_pgrid)then ! process U,V
      nvar3dh=6
   else
      nvar3dh=4
   endif
endif

if (make_micro) then
   nvar3dh=nvar3dh+8
endif

print*,'hinterp grid allocation 3d/3dpts = ',nvar3dh,lx*ly*nz*nvar3dh    
print*,'ctmax (predicted) = ',nz*nvar3dh    
print*,'nvar3dh / nvar3d = ',nvar3dh,nvar3d

allocate(hgrid(lx,ly,nvar3dh*nz))                                                

hgrid=rmsg

ct=1

hpsig    =>hgrid(1:lx,1:ly,ct:ct+nz-1); ct=ct+nz

hzsig    =>hgrid(1:lx,1:ly,ct:ct+nz-1); ct=ct+nz

htsig    =>hgrid(1:lx,1:ly,ct:ct+nz-1); ct=ct+nz
hmrsig   =>hgrid(1:lx,1:ly,ct:ct+nz-1); ct=ct+nz

if(.not. large_pgrid)then
   husig    =>hgrid(1:lx,1:ly,ct:ct+nz-1); ct=ct+nz
   hvsig    =>hgrid(1:lx,1:ly,ct:ct+nz-1); ct=ct+nz
   if(.not. large_ngrid)then   
      hwsig    =>hgrid(1:lx,1:ly,ct:ct+nz-1); ct=ct+nz
   endif
!  htkesig  =>hgrid(1:lx,1:ly,ct:ct+nz-1); ct=ct+nz
endif

if (make_micro) then
   hcldliqmr_sig =>hgrid(1:lx,1:ly,ct:ct+nz-1); ct=ct+nz
   hcldicemr_sig =>hgrid(1:lx,1:ly,ct:ct+nz-1); ct=ct+nz
   hrainmr_sig   =>hgrid(1:lx,1:ly,ct:ct+nz-1); ct=ct+nz
   hsnowmr_sig   =>hgrid(1:lx,1:ly,ct:ct+nz-1); ct=ct+nz
   hgraupelmr_sig=>hgrid(1:lx,1:ly,ct:ct+nz-1); ct=ct+nz
   hrefl_sig     =>hgrid(1:lx,1:ly,ct:ct+nz-1); ct=ct+nz
   hzdr_sig     =>hgrid(1:lx,1:ly,ct:ct+nz-1); ct=ct+nz
   hldr_sig     =>hgrid(1:lx,1:ly,ct:ct+nz-1); ct=ct+nz
endif

print*,'ctmax (actual) = ',ct - 1            
if(nz*nvar3dh .ne. ct-1)then
    write(6,*)' ERROR: ctmax actual is different from predicted'
    stop
endif

return
end subroutine

!===============================================================================

subroutine alloc_isobaric_grid

implicit none

integer :: ct

if(.not. large_pgrid)then ! state variables
   nvar3dout=8
else
   nvar3dout=0
endif

if (make_micro) then
   if(.not. large_pgrid)then
      nvar3dout=nvar3dout+9
   else
      nvar3dout=nvar3dout+1
   endif
endif

allocate(pgrid(lx,ly,nvar3dout*lz),name3d(nvar3dout*lz),units3d(nvar3dout*lz)  &
        ,lvltype3d(nvar3dout*lz),com3d(nvar3dout*lz),lvls3d(nvar3dout*lz))

if (out_grib) then
   allocate(gribitua(nvar3dout),paramua(nvar3dout),scalep10ua(nvar3dout))
endif

pgrid=rmsg

ct=1

if(.not. large_pgrid)then
   zprs  =>pgrid(1:lx,1:ly,ct:ct+lz-1); name3d(ct:ct+lz-1)='HT '; com3d(ct:ct+lz-1)='Geopotential Height'       ; lvls3d(ct:ct+lz-1)=nint(lprs); ct=ct+lz
   rhprs =>pgrid(1:lx,1:ly,ct:ct+lz-1); name3d(ct:ct+lz-1)='RH3'; com3d(ct:ct+lz-1)='Relative Humidity'         ; lvls3d(ct:ct+lz-1)=nint(lprs); ct=ct+lz
   tprs  =>pgrid(1:lx,1:ly,ct:ct+lz-1); name3d(ct:ct+lz-1)='T3 '; com3d(ct:ct+lz-1)='Temperature'               ; lvls3d(ct:ct+lz-1)=nint(lprs); ct=ct+lz
   shprs =>pgrid(1:lx,1:ly,ct:ct+lz-1); name3d(ct:ct+lz-1)='SH '; com3d(ct:ct+lz-1)='Specific Humidity'         ; lvls3d(ct:ct+lz-1)=nint(lprs); ct=ct+lz
   uprs  =>pgrid(1:lx,1:ly,ct:ct+lz-1); name3d(ct:ct+lz-1)='U3 '; com3d(ct:ct+lz-1)='U-component Wind'          ; lvls3d(ct:ct+lz-1)=nint(lprs); ct=ct+lz
   vprs  =>pgrid(1:lx,1:ly,ct:ct+lz-1); name3d(ct:ct+lz-1)='V3 '; com3d(ct:ct+lz-1)='V-component Wind'          ; lvls3d(ct:ct+lz-1)=nint(lprs); ct=ct+lz
!  tkeprs=>pgrid(1:lx,1:ly,ct:ct+lz-1); name3d(ct:ct+lz-1)='TKE'; com3d(ct:ct+lz-1)='Turbulent Kinetic Energy'  ; lvls3d(ct:ct+lz-1)=nint(lprs); ct=ct+lz
   wprs  =>pgrid(1:lx,1:ly,ct:ct+lz-1); name3d(ct:ct+lz-1)='W3 '; com3d(ct:ct+lz-1)='Vertical Velocity'         ; lvls3d(ct:ct+lz-1)=nint(lprs); ct=ct+lz
   omprs =>pgrid(1:lx,1:ly,ct:ct+lz-1); name3d(ct:ct+lz-1)='OM '; com3d(ct:ct+lz-1)='Pressure Vertical Velocity'; lvls3d(ct:ct+lz-1)=nint(lprs); ct=ct+lz
endif

if (make_micro) then
   if(.not. large_pgrid)then
      cldliqmr_prs =>pgrid(1:lx,1:ly,ct:ct+lz-1); name3d(ct:ct+lz-1)='LWC'; com3d(ct:ct+lz-1)='Cloud Liq.'     ; lvls3d(ct:ct+lz-1)=nint(lprs); ct=ct+lz
      cldicemr_prs =>pgrid(1:lx,1:ly,ct:ct+lz-1); name3d(ct:ct+lz-1)='ICE'; com3d(ct:ct+lz-1)='Cloud Ice'      ; lvls3d(ct:ct+lz-1)=nint(lprs); ct=ct+lz
      rainmr_prs   =>pgrid(1:lx,1:ly,ct:ct+lz-1); name3d(ct:ct+lz-1)='RAI'; com3d(ct:ct+lz-1)='Rain Conc.'     ; lvls3d(ct:ct+lz-1)=nint(lprs); ct=ct+lz
      snowmr_prs   =>pgrid(1:lx,1:ly,ct:ct+lz-1); name3d(ct:ct+lz-1)='SNO'; com3d(ct:ct+lz-1)='Snow Conc.'     ; lvls3d(ct:ct+lz-1)=nint(lprs); ct=ct+lz
      graupelmr_prs=>pgrid(1:lx,1:ly,ct:ct+lz-1); name3d(ct:ct+lz-1)='PIC'; com3d(ct:ct+lz-1)='Graupel Conc.'  ; lvls3d(ct:ct+lz-1)=nint(lprs); ct=ct+lz
      refl_prs     =>pgrid(1:lx,1:ly,ct:ct+lz-1); name3d(ct:ct+lz-1)='REF'; com3d(ct:ct+lz-1)='Radar Ref.'     ; lvls3d(ct:ct+lz-1)=nint(lprs); ct=ct+lz
      zdr_prs     =>pgrid(1:lx,1:ly,ct:ct+lz-1); name3d(ct:ct+lz-1)='ZDR'; com3d(ct:ct+lz-1)='Radar ZDR'       ; lvls3d(ct:ct+lz-1)=nint(lprs); ct=ct+lz
      ldr_prs     =>pgrid(1:lx,1:ly,ct:ct+lz-1); name3d(ct:ct+lz-1)='LDR'; com3d(ct:ct+lz-1)='Radar LDR'       ; lvls3d(ct:ct+lz-1)=nint(lprs); ct=ct+lz
      pcptype_prs  =>pgrid(1:lx,1:ly,ct:ct+lz-1); name3d(ct:ct+lz-1)='PTY'; com3d(ct:ct+lz-1)='Precip. Type'   ; lvls3d(ct:ct+lz-1)=nint(lprs); ct=ct+lz
   else ! allocate just reflectivity 
      refl_prs     =>pgrid(1:lx,1:ly,ct:ct+lz-1); name3d(ct:ct+lz-1)='REF'; com3d(ct:ct+lz-1)='Radar Ref.'     ; lvls3d(ct:ct+lz-1)=nint(lprs); ct=ct+lz
   endif
endif

if (out_grib) then
   ct=1
   gribitua(ct)=.true.; paramua(ct)=7  ; scalep10(ct)=0 ; ct=ct+1  ! zprs
   gribitua(ct)=.false.;                                ; ct=ct+1  ! rhprs
   gribitua(ct)=.true.; paramua(ct)=11 ; scalep10(ct)=2 ; ct=ct+1  ! tprs
   gribitua(ct)=.true.; paramua(ct)=51 ; scalep10(ct)=8 ; ct=ct+1  ! shprs
   gribitua(ct)=.true.; paramua(ct)=33 ; scalep10(ct)=1 ; ct=ct+1  ! uprs
   gribitua(ct)=.true.; paramua(ct)=34 ; scalep10(ct)=1 ; ct=ct+1  ! vprs
!   gribitua(ct)=.true.; paramua(ct)=158; scalep10(ct)=3 ; ct=ct+1  ! tkeprs
   gribitua(ct)=.true.; paramua(ct)=40 ; scalep10(ct)=3 ; ct=ct+1  ! wprs
   gribitua(ct)=.true.; paramua(ct)=39 ; scalep10(ct)=3 ; ct=ct+1  ! omprs
   if (make_micro) then
      gribitua(ct)=.true.; paramua(ct)=153; scalep10(ct)=6 ; ct=ct+1  ! cldliqmr_prs
      gribitua(ct)=.true.; paramua(ct)=178; scalep10(ct)=6 ; ct=ct+1  ! cldicemr_prs
      gribitua(ct)=.true.; paramua(ct)=170; scalep10(ct)=6 ; ct=ct+1  ! rainmr_prs
      gribitua(ct)=.true.; paramua(ct)=171; scalep10(ct)=6 ; ct=ct+1  ! snowmr_prs
      gribitua(ct)=.true.; paramua(ct)=179; scalep10(ct)=6 ; ct=ct+1  ! graupelmr_prs
      gribitua(ct)=.true.; paramua(ct)=128; scalep10(ct)=0 ; ct=ct+1  ! refl_prs
!      gribitua(ct)=.true.; paramua(ct)=128; scalep10(ct)=0 ; ct=ct+1  ! zdr_prs
!      gribitua(ct)=.true.; paramua(ct)=128; scalep10(ct)=0 ; ct=ct+1  ! ldr_prs
      gribitua(ct)=.true.; paramua(ct)=136; scalep10(ct)=0 ; ct=ct+1  ! pcptype_prs
   endif
endif

return
end subroutine

!===============================================================================

subroutine alloc_surface_grid

implicit none

integer :: ct

if (trim(mtype) /= 'st4') then
  nvar2dout=52
  if (make_micro) nvar2dout=nvar2dout+5
  if (make_firewx) nvar2dout=nvar2dout+7

  allocate(sgrid(lx,ly,nvar2dout),name2d(nvar2dout),units2d(nvar2dout)  &
        ,lvltype2d(nvar2dout),com2d(nvar2dout),lvls2d(nvar2dout))     

  if (out_grib) then
     allocate(gribit(nvar2dout),param(nvar2dout),leveltype(nvar2dout)  &
           ,level1(nvar2dout),level2(nvar2dout),timerange(nvar2dout),scalep10(nvar2dout))
  endif

  sgrid=rmsg

  ct=1

! The order of the first set of variables needs to match the order of 
!   the surface variables in the native grid.

  zsfc       =>sgrid(1:lx,1:ly,ct); name2d(ct)='TER'; com2d(ct)='Model Terrain'              ; ct=ct+1
  psfc       =>sgrid(1:lx,1:ly,ct); name2d(ct)='PSF'; com2d(ct)='Surface Pressure'           ; ct=ct+1
  tsfc       =>sgrid(1:lx,1:ly,ct); name2d(ct)='TSF'; com2d(ct)='Sfc Temperature'            ; ct=ct+1
  mrsfc      =>sgrid(1:lx,1:ly,ct); name2d(ct)='MSF'; com2d(ct)='Sfc Mixing Ratio'           ; ct=ct+1
  usfc       =>sgrid(1:lx,1:ly,ct); name2d(ct)='USF'; com2d(ct)='Sfc U-Component Wind'       ; ct=ct+1
  vsfc       =>sgrid(1:lx,1:ly,ct); name2d(ct)='VSF'; com2d(ct)='Sfc V-Component Wind'       ; ct=ct+1
  wsfc       =>sgrid(1:lx,1:ly,ct); name2d(ct)='WSF'; com2d(ct)='Sfc Vertical Velocity'      ; ct=ct+1
  ground_t   =>sgrid(1:lx,1:ly,ct); name2d(ct)='TGD'; com2d(ct)='Ground Temperature'         ; ct=ct+1
  pblhgt     =>sgrid(1:lx,1:ly,ct); name2d(ct)='BLH'; com2d(ct)='Boundary Layer Depth'       ; ct=ct+1
  pcp_tot    =>sgrid(1:lx,1:ly,ct); name2d(ct)='RTO'; com2d(ct)='Run-total Liq. Precip Accum'; ct=ct+1
  lwout      =>sgrid(1:lx,1:ly,ct); name2d(ct)='LWO'; com2d(ct)='Outgoing LW Radiation'      ; ct=ct+1
  swout      =>sgrid(1:lx,1:ly,ct); name2d(ct)='SWO'; com2d(ct)='Outgoing SW Radiation'      ; ct=ct+1
  lwdown     =>sgrid(1:lx,1:ly,ct); name2d(ct)='LWI'; com2d(ct)='Incoming LW Radiation'      ; ct=ct+1
  swdown     =>sgrid(1:lx,1:ly,ct); name2d(ct)='SWI'; com2d(ct)='Incoming SW Radiation'      ; ct=ct+1
  shflux     =>sgrid(1:lx,1:ly,ct); name2d(ct)='SHF'; com2d(ct)='Sensible Heat Flux'         ; ct=ct+1
  lhflux     =>sgrid(1:lx,1:ly,ct); name2d(ct)='LHF'; com2d(ct)='Latent Heat Flux'           ; ct=ct+1
  upflux     =>sgrid(1:lx,1:ly,ct); name2d(ct)='UMF'; com2d(ct)='Upslope Moisture Flux'      ; ct=ct+1
  bt11u      =>sgrid(1:lx,1:ly,ct); name2d(ct)='S8A'; com2d(ct)='11u Brightness Temperature' ; ct=ct+1
  simvis     =>sgrid(1:lx,1:ly,ct); name2d(ct)='SMV'; com2d(ct)='Albedo (Vis Satellite)'     ; ct=ct+1


  thetasfc   =>sgrid(1:lx,1:ly,ct); name2d(ct)='TH '; com2d(ct)='Sfc Potential Temperature'       ; ct=ct+1
  thetaesfc  =>sgrid(1:lx,1:ly,ct); name2d(ct)='THE'; com2d(ct)='Sfc Equiv. Potential Temperature'; ct=ct+1
  rhsfc      =>sgrid(1:lx,1:ly,ct); name2d(ct)='RH '; com2d(ct)='Sfc Relative Humidity'           ; ct=ct+1
  tdsfc      =>sgrid(1:lx,1:ly,ct); name2d(ct)='DSF'; com2d(ct)='Sfc Dewpoint Temperature'        ; ct=ct+1
  redp       =>sgrid(1:lx,1:ly,ct); name2d(ct)='P  '; com2d(ct)='Reduced Pressure'                ; ct=ct+1
  pmsl       =>sgrid(1:lx,1:ly,ct); name2d(ct)='SLP'; com2d(ct)='Sea-level Pressure'              ; ct=ct+1
  ztw0       =>sgrid(1:lx,1:ly,ct); name2d(ct)='TW0'; com2d(ct)='Height of wet-bulb = zero'       ; ct=ct+1
  ztw1       =>sgrid(1:lx,1:ly,ct); name2d(ct)='TW1'; com2d(ct)='Height of wet-bulb = 1.3'        ; ct=ct+1
  cldbase    =>sgrid(1:lx,1:ly,ct); name2d(ct)='LCB'; com2d(ct)='Cloud Base ASL'                  ; ct=ct+1
  cldtop     =>sgrid(1:lx,1:ly,ct); name2d(ct)='LCT'; com2d(ct)='Cloud Top ASL'                   ; ct=ct+1
  cldamt     =>sgrid(1:lx,1:ly,ct); name2d(ct)='LCV'; com2d(ct)='Cloud Opacity'                   ; ct=ct+1
  cldalb     =>sgrid(1:lx,1:ly,ct); name2d(ct)='CLA'; com2d(ct)='Cloud Albedo'                    ; ct=ct+1
  ceiling    =>sgrid(1:lx,1:ly,ct); name2d(ct)='CCE'; com2d(ct)='Cloud Ceiling AGL'               ; ct=ct+1
  intliqwater=>sgrid(1:lx,1:ly,ct); name2d(ct)='LIL'; com2d(ct)='Integrated Cloud Liquid'         ; ct=ct+1
  intcldice  =>sgrid(1:lx,1:ly,ct); name2d(ct)='LIC'; com2d(ct)='Integrated Cloud Ice'            ; ct=ct+1
! intgraupel =>sgrid(1:lx,1:ly,ct); name2d(ct)='LIG'; com2d(ct)='Integrated Graupel'              ; ct=ct+1
  totpcpwater=>sgrid(1:lx,1:ly,ct); name2d(ct)='TPW'; com2d(ct)='Total Precipitable Water'        ; ct=ct+1
  max_refl   =>sgrid(1:lx,1:ly,ct); name2d(ct)='LMR'; com2d(ct)='Composite Reflectivity'          ; ct=ct+1
  echo_tops  =>sgrid(1:lx,1:ly,ct); name2d(ct)='LMT'; com2d(ct)='Radar Echo Tops'                 ; ct=ct+1
  refl_sfc   =>sgrid(1:lx,1:ly,ct); name2d(ct)='LLR'; com2d(ct)='1km AGL Reflectivity'            ; ct=ct+1
  pcptype_sfc=>sgrid(1:lx,1:ly,ct); name2d(ct)='SPT'; com2d(ct)='Sfc Precip. Type'                ; ct=ct+1
  pcp_inc    =>sgrid(1:lx,1:ly,ct); name2d(ct)='R01'; com2d(ct)='Incremental Tot. Liq. Precip'    ; ct=ct+1
  snow_inc   =>sgrid(1:lx,1:ly,ct); name2d(ct)='S01'; com2d(ct)='Incremental Snow Depth'          ; ct=ct+1
  snow_tot   =>sgrid(1:lx,1:ly,ct); name2d(ct)='STO'; com2d(ct)='Run-total Snow Accum'            ; ct=ct+1
  srhel      =>sgrid(1:lx,1:ly,ct); name2d(ct)='LHE'; com2d(ct)='Storm Relative Helicity'         ; ct=ct+1
  uhel       =>sgrid(1:lx,1:ly,ct); name2d(ct)='UHE'; com2d(ct)='Updraft Helicity'                ; ct=ct+1
  cape       =>sgrid(1:lx,1:ly,ct); name2d(ct)='PBE'; com2d(ct)='CAPE'                            ; ct=ct+1
  cin        =>sgrid(1:lx,1:ly,ct); name2d(ct)='NBE'; com2d(ct)='CIN'                             ; ct=ct+1
  liftedind  =>sgrid(1:lx,1:ly,ct); name2d(ct)='LI '; com2d(ct)='Lifted Index'                    ; ct=ct+1
  visibility =>sgrid(1:lx,1:ly,ct); name2d(ct)='VIS'; com2d(ct)='Sfc. Visibility'                 ; ct=ct+1
  heatind    =>sgrid(1:lx,1:ly,ct); name2d(ct)='HI '; com2d(ct)='Heat Index'                      ; ct=ct+1
  u80        =>sgrid(1:lx,1:ly,ct); name2d(ct)='U80'; com2d(ct)='U-component Wind at 80m'         ; ct=ct+1
  v80        =>sgrid(1:lx,1:ly,ct); name2d(ct)='V80'; com2d(ct)='V-component Wind at 80m'         ; ct=ct+1

  if (make_micro) then
     clwmrsfc  =>sgrid(1:lx,1:ly,ct); name2d(ct)='SCL'; com2d(ct)='Sfc Cloud Liq Water MR'; ct=ct+1
     icemrsfc  =>sgrid(1:lx,1:ly,ct); name2d(ct)='SIC'; com2d(ct)='Sfc Cloud Ice MR'      ; ct=ct+1
     snowmrsfc =>sgrid(1:lx,1:ly,ct); name2d(ct)='SSN'; com2d(ct)='Sfc Prec. Snow MR'     ; ct=ct+1
     rainmrsfc =>sgrid(1:lx,1:ly,ct); name2d(ct)='SRN'; com2d(ct)='Sfc Prec. Rain MR'     ; ct=ct+1
     graupmrsfc=>sgrid(1:lx,1:ly,ct); name2d(ct)='SGR'; com2d(ct)='Sfc Prec. Graupel MR'  ; ct=ct+1
  endif

  if (make_firewx) then
     upbl     =>sgrid(1:lx,1:ly,ct); name2d(ct)='UPB'; com2d(ct)='U-component Wind in PBL'; ct=ct+1
     vpbl     =>sgrid(1:lx,1:ly,ct); name2d(ct)='VPB'; com2d(ct)='V-component Wind in PBL'; ct=ct+1
     vnt_index=>sgrid(1:lx,1:ly,ct); name2d(ct)='VNT'; com2d(ct)='Ventilation Index'      ; ct=ct+1
     ham_index=>sgrid(1:lx,1:ly,ct); name2d(ct)='HAM'; com2d(ct)='Mid-Level Haines Index' ; ct=ct+1
     hah_index=>sgrid(1:lx,1:ly,ct); name2d(ct)='HAH'; com2d(ct)='High-Level Haines Index'; ct=ct+1
     fwi_index=>sgrid(1:lx,1:ly,ct); name2d(ct)='FWI'; com2d(ct)='Fosberg Fire Wx Index'  ; ct=ct+1
     fwx_index=>sgrid(1:lx,1:ly,ct); name2d(ct)='FWX'; com2d(ct)='LAPS/Kelsch Fire Wx Index'  ; ct=ct+1
  endif

  if (out_grib) then
     ct=1
     gribit(ct)=.true.; param(ct)=7  ; leveltype(ct)=1  ; level1(ct)=0 ; level2(ct)=0; timerange(ct)=0; scalep10(ct)=0 ; ct=ct+1  ! zsfc
     gribit(ct)=.true.; param(ct)=1  ; leveltype(ct)=1  ; level1(ct)=0 ; level2(ct)=0; timerange(ct)=0; scalep10(ct)=1 ; ct=ct+1  ! psfc
     gribit(ct)=.true.; param(ct)=11 ; leveltype(ct)=105; level1(ct)=2 ; level2(ct)=0; timerange(ct)=0; scalep10(ct)=1 ; ct=ct+1  ! tsfc
     gribit(ct)=.false.                                                                                                ; ct=ct+1  ! mrsfc
     gribit(ct)=.true.; param(ct)=33 ; leveltype(ct)=105; level1(ct)=10; level2(ct)=0; timerange(ct)=0; scalep10(ct)=1 ; ct=ct+1  ! usfc
     gribit(ct)=.true.; param(ct)=34 ; leveltype(ct)=105; level1(ct)=10; level2(ct)=0; timerange(ct)=0; scalep10(ct)=1 ; ct=ct+1  ! vsfc
     gribit(ct)=.true.; param(ct)=40 ; leveltype(ct)=105; level1(ct)=10; level2(ct)=0; timerange(ct)=0; scalep10(ct)=4 ; ct=ct+1  ! wsfc
     gribit(ct)=.true.; param(ct)=11 ; leveltype(ct)=1  ; level1(ct)=0 ; level2(ct)=0; timerange(ct)=0; scalep10(ct)=1 ; ct=ct+1  ! ground_t
     gribit(ct)=.true.; param(ct)=221; leveltype(ct)=1  ; level1(ct)=0 ; level2(ct)=0; timerange(ct)=0; scalep10(ct)=0 ; ct=ct+1  ! pblhgt
     gribit(ct)=.true.; param(ct)=142; leveltype(ct)=1  ; level1(ct)=0 ; level2(ct)=0; timerange(ct)=4; scalep10(ct)=1 ; ct=ct+1  ! pcp_tot
     gribit(ct)=.true.; param(ct)=212; leveltype(ct)=8  ; level1(ct)=0 ; level2(ct)=0; timerange(ct)=0; scalep10(ct)=0 ; ct=ct+1  ! lwout
     gribit(ct)=.true.; param(ct)=211; leveltype(ct)=8  ; level1(ct)=0 ; level2(ct)=0; timerange(ct)=0; scalep10(ct)=0 ; ct=ct+1  ! swout
     gribit(ct)=.true.; param(ct)=112; leveltype(ct)=1  ; level1(ct)=0 ; level2(ct)=0; timerange(ct)=0; scalep10(ct)=0 ; ct=ct+1  ! lwdown
     gribit(ct)=.true.; param(ct)=111; leveltype(ct)=1  ; level1(ct)=0 ; level2(ct)=0; timerange(ct)=0; scalep10(ct)=0 ; ct=ct+1  ! swdown
     gribit(ct)=.true.; param(ct)=122; leveltype(ct)=1  ; level1(ct)=0 ; level2(ct)=0; timerange(ct)=0; scalep10(ct)=0 ; ct=ct+1  ! shflux
     gribit(ct)=.true.; param(ct)=121; leveltype(ct)=1  ; level1(ct)=0 ; level2(ct)=0; timerange(ct)=0; scalep10(ct)=0 ; ct=ct+1  ! lhflux
!     gribit(ct)=.true.; param(ct)=121; leveltype(ct)=1  ; level1(ct)=0 ; level2(ct)=0; timerange(ct)=0; scalep10(ct)=0 ; ct=ct+1  ! upflux

     gribit(ct)=.true.; param(ct)=13 ; leveltype(ct)=105; level1(ct)=2 ; level2(ct)=0; timerange(ct)=0; scalep10(ct)=1 ; ct=ct+1  ! thetasfc
     gribit(ct)=.true.; param(ct)=14 ; leveltype(ct)=105; level1(ct)=2 ; level2(ct)=0; timerange(ct)=0; scalep10(ct)=1 ; ct=ct+1  ! thetaesfc
     gribit(ct)=.true.; param(ct)=52 ; leveltype(ct)=105; level1(ct)=2 ; level2(ct)=0; timerange(ct)=0; scalep10(ct)=1 ; ct=ct+1  ! rhsfc
     gribit(ct)=.true.; param(ct)=17 ; leveltype(ct)=105; level1(ct)=2 ; level2(ct)=0; timerange(ct)=0; scalep10(ct)=1 ; ct=ct+1  ! tdsfc
     gribit(ct)=.false.                                                                                                ; ct=ct+1  ! redp
     gribit(ct)=.true.; param(ct)=2  ; leveltype(ct)=102; level1(ct)=0 ; level2(ct)=0; timerange(ct)=0; scalep10(ct)=1 ; ct=ct+1  ! pmsl
     gribit(ct)=.false.                                                                                                ; ct=ct+1  ! ztw0
     gribit(ct)=.false.                                                                                                ; ct=ct+1  ! ztw1
     gribit(ct)=.true.; param(ct)=138; leveltype(ct)=102; level1(ct)=0 ; level2(ct)=0; timerange(ct)=0; scalep10(ct)=0 ; ct=ct+1  ! cldbase
     gribit(ct)=.true.; param(ct)=139; leveltype(ct)=102; level1(ct)=0 ; level2(ct)=0; timerange(ct)=0; scalep10(ct)=0 ; ct=ct+1  ! cldtop
     gribit(ct)=.true.; param(ct)=71 ; leveltype(ct)=1  ; level1(ct)=0 ; level2(ct)=0; timerange(ct)=0; scalep10(ct)=0 ; ct=ct+1  ! cldamt
     gribit(ct)=.true.; param(ct)=137; leveltype(ct)=1  ; level1(ct)=0 ; level2(ct)=0; timerange(ct)=0; scalep10(ct)=0 ; ct=ct+1  ! ceiling
     gribit(ct)=.false.                                                                                                ; ct=ct+1  ! intliqwater
     gribit(ct)=.true.; param(ct)=54 ; leveltype(ct)=200; level1(ct)=0 ; level2(ct)=0; timerange(ct)=0; scalep10(ct)=0 ; ct=ct+1  ! totpcpwater
     gribit(ct)=.true.; param(ct)=129; leveltype(ct)=1  ; level1(ct)=0 ; level2(ct)=0; timerange(ct)=0; scalep10(ct)=0 ; ct=ct+1  ! max_refl
     gribit(ct)=.true.; param(ct)=130; leveltype(ct)=102; level1(ct)=0 ; level2(ct)=0; timerange(ct)=0; scalep10(ct)=0 ; ct=ct+1  ! echo_tops
     gribit(ct)=.true.; param(ct)=128; leveltype(ct)=1  ; level1(ct)=0 ; level2(ct)=0; timerange(ct)=0; scalep10(ct)=0 ; ct=ct+1  ! refl_sfc
     gribit(ct)=.true.; param(ct)=136; leveltype(ct)=1  ; level1(ct)=0 ; level2(ct)=0; timerange(ct)=0; scalep10(ct)=0 ; ct=ct+1  ! pcptype_sfc
     gribit(ct)=.true.; param(ct)=61 ; leveltype(ct)=1  ; level1(ct)=0 ; level2(ct)=0; timerange(ct)=4; scalep10(ct)=1 ; ct=ct+1  ! pcp_inc
     gribit(ct)=.true.; param(ct)=66 ; leveltype(ct)=1  ; level1(ct)=0 ; level2(ct)=0; timerange(ct)=4; scalep10(ct)=4 ; ct=ct+1  ! snow_inc
     gribit(ct)=.true.; param(ct)=141; leveltype(ct)=1  ; level1(ct)=0 ; level2(ct)=0; timerange(ct)=4; scalep10(ct)=4 ; ct=ct+1  ! snow_tot
     gribit(ct)=.true.; param(ct)=190; leveltype(ct)=1  ; level1(ct)=0 ; level2(ct)=0; timerange(ct)=0; scalep10(ct)=0 ; ct=ct+1  ! srhel
     gribit(ct)=.true.; param(ct)=157; leveltype(ct)=1  ; level1(ct)=0 ; level2(ct)=0; timerange(ct)=0; scalep10(ct)=0 ; ct=ct+1  ! cape
     gribit(ct)=.true.; param(ct)=156; leveltype(ct)=1  ; level1(ct)=0 ; level2(ct)=0; timerange(ct)=0; scalep10(ct)=0 ; ct=ct+1  ! cin
     gribit(ct)=.true.; param(ct)=131; leveltype(ct)=1  ; level1(ct)=0 ; level2(ct)=0; timerange(ct)=0; scalep10(ct)=1 ; ct=ct+1  ! liftedind
     gribit(ct)=.true.; param(ct)=20 ; leveltype(ct)=1  ; level1(ct)=0 ; level2(ct)=0; timerange(ct)=0; scalep10(ct)=-3; ct=ct+1  ! visibility
     gribit(ct)=.false.                                                                                                ; ct=ct+1  ! heatind

     if (make_micro) then
        gribit(ct)=.true.; param(ct)=153; leveltype(ct)=1  ; level1(ct)=0 ; level2(ct)=0; timerange(ct)=0; scalep10(ct)=4; ct=ct+1  ! clwmrsfc
        gribit(ct)=.true.; param(ct)=178; leveltype(ct)=1  ; level1(ct)=0 ; level2(ct)=0; timerange(ct)=0; scalep10(ct)=4; ct=ct+1  ! icemrsfc
        gribit(ct)=.true.; param(ct)=171; leveltype(ct)=1  ; level1(ct)=0 ; level2(ct)=0; timerange(ct)=0; scalep10(ct)=4; ct=ct+1  ! snowmrsfc
        gribit(ct)=.true.; param(ct)=170; leveltype(ct)=1  ; level1(ct)=0 ; level2(ct)=0; timerange(ct)=0; scalep10(ct)=4; ct=ct+1  ! rainmrsfc
        gribit(ct)=.true.; param(ct)=179; leveltype(ct)=1  ; level1(ct)=0 ; level2(ct)=0; timerange(ct)=0; scalep10(ct)=4; ct=ct+1  ! graupmrsfc
     endif

     if (make_firewx) then
        gribit(ct)=.false.                                                                                               ; ct=ct+1  ! upbl
        gribit(ct)=.false.                                                                                               ; ct=ct+1  ! vpbl
        gribit(ct)=.false.                                                                                               ; ct=ct+1  ! vnt_index
        gribit(ct)=.false.                                                                                               ; ct=ct+1  ! ham_index
        gribit(ct)=.false.                                                                                               ; ct=ct+1  ! hah_index
        gribit(ct)=.false.                                                                                               ; ct=ct+1  ! fwi_index
        gribit(ct)=.false.                                                                                               ; ct=ct+1  ! fwx_index
     endif
  endif
else
  nvar2dout=1

  allocate(sgrid(lx,ly,nvar2dout),name2d(nvar2dout),units2d(nvar2dout)  &
        ,lvltype2d(nvar2dout),com2d(nvar2dout),lvls2d(nvar2dout))     

  if (out_grib) then
     allocate(gribit(nvar2dout),param(nvar2dout),leveltype(nvar2dout)  &
           ,level1(nvar2dout),level2(nvar2dout),timerange(nvar2dout),scalep10(nvar2dout))
  endif

  sgrid=rmsg

  ct=1
  pcp_inc    =>sgrid(1:lx,1:ly,ct); name2d(ct)='R01'; com2d(ct)='Incremental Tot. Liq. Precip'    ; ct=ct+1
  if (out_grib) then
   ct=1
   gribit(ct)=.true.; param(ct)=61 ; leveltype(ct)=1  ; level1(ct)=0 ; level2(ct)=0; timerange(ct)=4; scalep10(ct)=1 ; ct=ct+1  ! pcp_inc
  endif
endif

return
end subroutine

!===============================================================================

subroutine dealloc_grid(gtype)

implicit none

character(len=*) :: gtype

select case(trim(gtype))
   case('native')
      deallocate(ngrid,nlat,nlon)
   case('horiz')
      deallocate(hgrid)
   case('isobaric')
      deallocate(pgrid,lprs,lprsl,name3d,units3d,lvltype3d,com3d,lvls3d)
   case('surface')
      deallocate(sgrid,llat,llon,name2d,units2d,lvltype2d,com2d,lvls2d)
      if (out_grib) deallocate(gribit,param,leveltype  &
                              ,level1,level2,timerange,scalep10)
end select

return
end subroutine

end module   
