subroutine lfm_namelist(lfmprd_dir) 

! Read namelist variables.
! Default values are set in lfmgrid

use lfmgrid

implicit none

integer :: istatus

character(len=*) :: lfmprd_dir
character(len=256) :: nlfile

namelist/lfmpost_nl/out_cdf,out_grib,out_v5d            &
                   ,make_micro,make_firewx,make_points  &
                   ,verbose,realtime,write_to_lapsdir   &
                   ,make_donefile                       &
                   ,precip_dt,c_m2z,advection_time      &
                   ,n3d_pts_thr,p3d_pts_thr             &
                   ,wrf_version
                   

nlfile=trim(lfmprd_dir)//'/../static/lfmpost.nl'
open(1,file=trim(nlfile),form='formatted',status='old',iostat=istatus)
if (istatus /= 0) then
   print*,'Error opening namelist file: ',trim(nlfile)
   stop
endif
read(1,nml=lfmpost_nl)
close(1)

return
end

!===============================================================================

subroutine get_native_dims(mtype,filename,nx,ny,nz,istatus)

use lfmgrid, only : n3d_pts_thr, large_ngrid 

implicit none

integer :: nx,ny,nz,istatus

character(len=*) :: mtype,filename

! istatus: 1=good return, 0=ERROR/bad return
istatus = 1  ! assume good return

select case(trim(mtype))
   case('mm5')
      call get_mm5_dims(filename,nx,ny,nz)
   case('wrf')
      call get_wrf_dims(filename,nx,ny,nz,istatus)
   case('nmm')
      call get_nmm_dims(filename,nx,ny,nz)
   case('st4')
      call get_st4_dims(filename,nx,ny,nz)
end select

! Set native grid size threshold based on an 8GB machine
write(6,*)' # of 3D native grid points, large_ngrid thresh = ',nx*ny*nz,n3d_pts_thr
if(nx*ny*nz > n3d_pts_thr)then
   write(6,*)' Large native grid - process reduced set of 3-D fields'
   large_ngrid = .true.
!  large_pgrid = .true.
endif

return
end

!===============================================================================

subroutine get_laps_static(laps_data_root)

use lfmgrid

implicit none

include 'netcdf.inc'

integer, parameter :: nprmax=150
integer :: icode,ncid,nid,nk_laps,istatus,i

real, dimension(nprmax) :: pressures

character(len=256) :: stlaps
character(len=132) :: gridtype
character(len=*) :: laps_data_root

namelist /pressures_nl/pressures

stlaps=trim(laps_data_root)//'/static/static.nest7grid'
icode=nf_open(trim(stlaps),nf_nowrite,ncid)
if (icode /= 0) then
   print*,'LAPS static file not found: ',trim(stlaps)
   stop
endif
icode=nf_inq_dimid(ncid,'x',nid)
icode=nf_inq_dimlen(ncid,nid,lx)
icode=nf_inq_dimid(ncid,'y',nid)
icode=nf_inq_dimlen(ncid,nid,ly)
icode=nf_inq_varid(ncid,'grid_type',nid)
icode=nf_get_var_text(ncid,nid,gridtype)
icode=nf_inq_varid(ncid,'grid_spacing',nid)
icode=nf_get_var_real(ncid,nid,grid_spacing)
icode=nf_inq_varid(ncid,'Latin1',nid)
icode=nf_get_var_real(ncid,nid,truelat1)
icode=nf_inq_varid(ncid,'Latin2',nid)
icode=nf_get_var_real(ncid,nid,truelat2)
icode=nf_inq_varid(ncid,'LoV',nid)
icode=nf_get_var_real(ncid,nid,stdlon)
if (stdlon > 180.) stdlon=stdlon-360.

allocate(llat(lx,ly),llon(lx,ly))
icode=nf_inq_varid(ncid,'lat',nid)
icode=nf_get_var_real(ncid,nid,llat)
icode=nf_inq_varid(ncid,'lon',nid)
icode=nf_get_var_real(ncid,nid,llon)
icode=nf_close(ncid)

if (gridtype(1:5) == 'polar') then
   projection='POLAR STEREOGRAPHIC'
elseif (gridtype(1:14) == 'secant lambert') then
   projection='LAMBERT CONFORMAL'
elseif (gridtype(1:18) == 'tangential lambert') then
   projection='LAMBERT CONFORMAL'
elseif (gridtype(1:18) == 'mercator') then
   projection='MERCATOR'
else
   print*,'Unrecognized LAPS grid type: ',trim(gridtype)
   stop
endif

pressures=rmsg
stlaps=trim(laps_data_root)//'/static/pressures.nl'
open(2,file=trim(stlaps),status='old',err=900)
read(2,pressures_nl,err=901)
close(2)

! Determine number of isobaric levels by checking pressure data.
!   (Assume there is at least one level, and that data is ordered correctly)

do lz=1,nprmax
   if (pressures(lz+1) > 200000.) exit
enddo

allocate(lprs(lz),lprsl(lz))
lprs(1:lz)=pressures(1:lz)
lprsl(1:lz)=alog(lprs(1:lz))

write(6,*)' # of 3D LAPS grid points, large_pgrid thresh = ',lx*ly*lz,p3d_pts_thr
if(lx*ly*lz > p3d_pts_thr)then
   write(6,*)' Large LAPS grid - process reduced set of 3-D fields'
   large_pgrid = .true.
!  large_ngrid = .true.
endif

return

900 continue
print*,'Could not open laps namelist file: ',trim(stlaps)
stop

901 continue
print*,'Error reading laps namelist file: ',trim(stlaps)
stop

end

!===============================================================================

subroutine set_laps_projection

use lfmgrid

implicit none

select case(trim(projection))
   case('LAMBERT CONFORMAL')
      call map_set(proj_lc,llat(1,1),llon(1,1),grid_spacing  &
                  ,stdlon,truelat1,truelat2,lx,ly,proj)
   case('POLAR STEREOGRAPHIC')
      call map_set(proj_ps,llat(1,1),llon(1,1),grid_spacing  &
                  ,stdlon,truelat1,truelat2,lx,ly,proj)
   case('MERCATOR')
      call map_set(proj_merc,llat(1,1),llon(1,1),grid_spacing  &
                  ,stdlon,truelat1,truelat2,lx,ly,proj)
end select

return
end

!===============================================================================

subroutine fill_native_grid

use lfmgrid

implicit none

! Fill native grids.

select case(trim(mtype))
   case('mm5')
      call fill_mm5_grid
   case('wrf')
      call fill_wrf_grid
   case('nmm')
      call fill_nmm_grid
   case('st4')
      call fill_st4_grid
end select

! Set map utilities for use by horizontal interpolation.

select case(trim(nprojection))
   case('LAMBERT CONFORMAL')
      call map_set(proj_lc,nlat(1,1),nlon(1,1),ngrid_spacingx  &
                  ,nstdlon,ntruelat1,ntruelat2,nx,ny,proj)
   case('POLAR STEREOGRAPHIC')
      call map_set(proj_ps,nlat(1,1),nlon(1,1),ngrid_spacingx  &
                  ,nstdlon,ntruelat1,ntruelat2,nx,ny,proj)
   case('MERCATOR')
      call map_set(proj_merc,nlat(1,1),nlon(1,1),ngrid_spacingx  &
                  ,nstdlon,ntruelat1,ntruelat2,nx,ny,proj)
   case('ROTATED LAT-LON')
end select

return
end

!===============================================================================

subroutine lfm_derived

use lfmgrid
use constants
use cloud_rad ! Cloud Radiation and Microphysics Parameters

implicit none

integer :: i,j,k,k500,k700,k850,k1000
real :: dz,dtdz,dewpt,relhum,potential_temp,eq_potential_temp,heatindex  &
       ,redp_lvl=1500.
real, allocatable, dimension(:,:) :: pcp_06,pcp_init,fallen_precip_type
real, allocatable, dimension(:,:,:) :: htdsig,hrhsig,hthetasig,hthetaesig,hzsigf  &
                                      ,pcpmr_sig,rhodrysig,tvsig,rhomoistsig

!beka
real :: dx(lx,ly),dy(lx,ly)  
real :: r_missing_data, aglhgt
real :: stefan_boltzmann, b_olr, eff_emissivity, b_tau
real :: bt_flux_equiv
real :: coeff
real :: a1,b1,c1,d1,e1,alpha(lx,ly)
real :: const_lwp,const_iwp,const_rwp,const_swp,const_gwp
real :: const_lwp_bks,const_iwp_bks,const_rwp_bks,const_swp_bks,const_gwp_bks
real :: clwc2alpha,cice2alpha,rain2alpha,snow2alpha,pice2alpha

!beka

character*256  directory
real           rspacing_dum
character*125   comment_2d
character*10    units_2d
character*50    ext
integer        len_dir,istatus
character*3    var_2d
real :: ldf(lx,ly),lat(lx,ly),lon(lx,ly),avg(lx,ly),static_albedo(lx,ly)
real :: windspeed(lx,ly),soil_moist(lx,ly),snow_cover(lx,ly)
real :: intrain(lx,ly),intsnow(lx,ly),intgraupel(lx,ly)
real :: const_pcpadj

integer ::        ismoist,isnow
integer ::        status

integer I4_elapsed, ishow_timer                    

!cj Variables to compute omega, added 6/21/2007
real :: tvprs

! Variables for low level reflectivity
integer :: klow, khigh
real :: ht_1km, frack
real :: umean(lx,ly),vmean(lx,ly),ustorm(lx,ly),vstorm(lx,ly)
real :: array_buf(lx,ly),array_out(lx,ly)

real :: ghi_ratio(lx,ly)
real :: arg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Beka!!!!!!!!!!!!! obtaining ldf, lat and lon !!!!!!!!!!!!!!!!!!!!

        ext = 'nest7grid'

        call get_directory(ext,directory,len_dir)
        var_2d='ldf'
        call rd_laps_static (directory,ext,lx,ly,1,var_2d, &
                             units_2d,comment_2d,ldf,      &
                             rspacing_dum,istatus)
        var_2d='avg'
        call rd_laps_static (directory,ext,lx,ly,1,var_2d, &
                             units_2d,comment_2d,avg,      &
                             rspacing_dum,istatus)

        var_2d='lat'
        call rd_laps_static (directory,ext,lx,ly,1,var_2d, &
                             units_2d,comment_2d,lat,      &
                             rspacing_dum,istatus)

        var_2d='lon'
        call rd_laps_static (directory,ext,lx,ly,1,var_2d, &
                             units_2d,comment_2d,lon,      &
                             rspacing_dum,istatus)

        call get_static_field_interp('albedo',laps_valtime,lx,ly &
                                    ,static_albedo,istatus)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Extrapolate a surface temp if not available.

if(.not. large_ngrid)then
  if (maxval(tsfc) > 1000.) then
   if (verbose) then
      print*,' '
      print*,'tsfc not available...will extrapolate from lowest level...'
   endif
   do j=1,ly
   do i=1,lx

! Compute dz between lowest sigma layer and 2m level

      dz=hzsig(i,j,1)-(zsfc(i,j)+2.)

! Compute the lapse rate of temp for the two lowest sigma levels.  
! Use of the hypsometric equation did not work well for these thin layers.

      dtdz=(htsig(i,j,2)-htsig(i,j,1))  &
          /(hzsig(i,j,2)-hzsig(i,j,1))
      tsfc(i,j)=htsig(i,j,1)-dtdz*dz
   enddo
   enddo
  endif
endif

! Fill other missing surface fields from lowest model level.

print *,'Min/Max mrsfc (model sfc grids) =',minval(mrsfc),maxval(mrsfc)
print *,'Substitute lowest sigma level for surface mixing ratio'

!if (maxval(mrsfc) > 1000.)mrsfc(:,:)=hmrsig(:,:,1)
if (.true.)               mrsfc(:,:)=hmrsig(:,:,1)

print *,'Min/Max mrsfc (model sfc grids) =',minval(mrsfc),maxval(mrsfc)

if(l_process_uv)then
  if (maxval(usfc) > 1000.) usfc(:,:)=husig(:,:,1)
  if (maxval(vsfc) > 1000.) vsfc(:,:)=hvsig(:,:,1)
endif
if(l_process_w)then
  if (maxval(wsfc) > 1000.) wsfc(:,:)=hwsig(:,:,1)
endif

! Other derived surface fields.
where(mrsfc < zero_thresh) mrsfc = zero_thresh

if (maxval(mrsfc) < 1000.)then
  do j=1,ly
  do i=1,lx
   rhsfc(i,j)=min(relhum(tsfc(i,j),mrsfc(i,j),psfc(i,j)),1.)
   tdsfc(i,j)=dewpt(tsfc(i,j),rhsfc(i,j))
   thetasfc(i,j)=potential_temp(tsfc(i,j),psfc(i,j))
   thetaesfc(i,j)=eq_potential_temp(tsfc(i,j),psfc(i,j),mrsfc(i,j),rhsfc(i,j))
   rhsfc(i,j)=rhsfc(i,j)*100.
  enddo
  enddo
endif

! Generate reduced pressure for LAPS usage.

if(.not. large_pgrid)then
 if (verbose) then
   print*,' '
   print*,'Reducing pressure to ',redp_lvl, ' meters'
 endif
 call interp_press_to_z(lprs,zprs,redp_lvl,redp,lx,ly,lz)

! Use same routine to interpolate sea-level pressure.  This method
! is used in lieu of original reduction routine, because it will
! keep the MSL field consistent with the height field, which has been
! properly reduced.  It also produces a bit smoother field over the mountains.

 if (verbose) print*,'Reducing pressure to MSL...'
 call interp_press_to_z(lprs,zprs,0.,pmsl,lx,ly,lz)
endif

! Microphysical surface variables.

if (make_micro) then
 where(hmrsig < zero_thresh) hmrsig=zero_thresh
 where(hcldliqmr_sig < zero_thresh) hcldliqmr_sig=0.0
 where(hcldicemr_sig < zero_thresh) hcldicemr_sig=0.0
 where(hsnowmr_sig < zero_thresh) hsnowmr_sig=0.0
 where(hrainmr_sig < zero_thresh) hrainmr_sig=0.0
 where(hgraupelmr_sig < zero_thresh) hgraupelmr_sig=0.0
 clwmrsfc(:,:)=hcldliqmr_sig(:,:,1)
 icemrsfc(:,:)=hcldicemr_sig(:,:,1)
 snowmrsfc(:,:)=hsnowmr_sig(:,:,1)
 rainmrsfc(:,:)=hrainmr_sig(:,:,1)
 graupmrsfc(:,:)=hgraupelmr_sig(:,:,1)

!if(.not. large_ngrid)then
 if(.true.)then                  

! Generate cloud fields.

   cldbase=rmsg
   cldtop=rmsg
   ceiling=rmsg  ! unlimited ceiling
   if (maxval(hcldliqmr_sig) < rmsg .and. maxval(hcldicemr_sig) < rmsg .and.  &
       maxval(hsnowmr_sig) < rmsg) then
      cldamt=0.
      call lfmclouds(lx,ly,nz,ngrid_spacingx,hcldliqmr_sig,hcldicemr_sig  &
                    ,hsnowmr_sig,hzsig,zsfc,cldbase,cldtop,ceiling,cldamt)
   else
      cldamt=rmsg
   endif

! Generate integrated hydrometeor fields

   allocate(pcpmr_sig(lx,ly,nz),rhodrysig(lx,ly,nz)) ! for possible use
   pcpmr_sig=0.
   if (maxval(hsnowmr_sig) < rmsg) pcpmr_sig=pcpmr_sig+hsnowmr_sig
   if (maxval(hrainmr_sig) < rmsg) pcpmr_sig=pcpmr_sig+hrainmr_sig
   if (maxval(hgraupelmr_sig) < rmsg) pcpmr_sig=pcpmr_sig+hgraupelmr_sig
   where(pcpmr_sig < zero_thresh) pcpmr_sig=0.
   rhodrysig=hpsig/(r*htsig)
   call lfm_integrated_liquid(lx,ly,nz,hcldliqmr_sig,hcldicemr_sig,hmrsig,&
       hrainmr_sig,hsnowmr_sig,hgraupelmr_sig,rhodrysig,hzsig,zsfc,  &
       intliqwater,intcldice,totpcpwater,intrain,intsnow,intgraupel,rmsg)
 endif ! large_ngrid

 if(.true.)then ! supercedes binary cldamt calculated earlier
                ! units of integrated vertical quantities (e.g. intliqwater) is meters
                ! for liquid we can augment this noting that tau = (3 * intliqwater * rholiq) / (2 * rho * radius)
                ! cldamt should be related to opacity and optical depth
                ! we might use a new module for this derived from these files
                ! in  'src/lib/cloud/get_cloud_rad.f90':
                ! and 'src/lib/modules/module_cloud_rad.f90'

  const_lwp = (1.5 * rholiq) / (rholiq     * reff_clwc)
  const_lwp_bks = const_lwp * bksct_eff_clwc

  const_iwp = (1.5 * rholiq) / (rholiq     * reff_cice)
  const_iwp_bks = const_iwp * bksct_eff_cice

  const_rwp = (1.5 * rholiq) / (rholiq     * reff_rain)
  const_rwp_bks = const_rwp * bksct_eff_rain

  const_swp = (1.5 * rholiq) / (rhosnow    * reff_snow)
  const_swp_bks = const_swp * bksct_eff_snow

  const_gwp = (1.5 * rholiq) / (rhograupel * reff_graupel)
  const_gwp_bks = const_gwp * bksct_eff_graupel

  write(6,*)' simvis constants ',const_lwp_bks,const_iwp_bks,const_rwp_bks,const_swp_bks,const_gwp_bks

  do j=1,ly
    do i=1,lx
      if(intliqwater(i,j) < rmsg .AND. intcldice(i,j) < rmsg)then
!        cldamt(i,j) = (intliqwater(i,j) * 100.) ** 0.3333333 ! this might work empirically if small integrated liq/ice
!        cldamt(i,j) = min(cldamt(i,j),1.0)                   ! values have smaller mean diameters, so the dependence
!                                                             ! is a weaker one
!        liquid water is assumed to have smaller particles contributing to a larger constant
!        Cloud amount is opacity of cloud liquid and cloud ice hydrometeors
         cldamt(i,j) = 1. - (exp( -(const_lwp * intliqwater(i,j) + const_iwp * intcldice(i,j))) ) 

!        Rain, Snow, and Graupel are added for the cldalb & simvis computation
         b_tau = const_lwp_bks * intliqwater(i,j) + const_iwp_bks * intcldice(i,j) + &
                 const_rwp_bks * intrain(i,j)     + const_swp_bks * intsnow(i,j)   + const_gwp_bks * intgraupel(i,j) 
         cldalb(i,j) = b_tau / (b_tau + 1.)
         simvis(i,j) = cldalb(i,j) + (1.-cldalb(i,j))**2 * (static_albedo(i,j)/(1.-cldalb(i,j)*static_albedo(i,j)))
      endif
    enddo ! i
  enddo ! j
 endif ! true

 if(.not. large_ngrid)then
   if (verbose) then
      print*,' '
      print *,'Called integrated_liquid.'
      print *,'Min/Max pcpmr_sig =',minval(pcpmr_sig),maxval(pcpmr_sig)
      print *,'Min/Max mrsig =',minval(hmrsig),maxval(hmrsig)
      print *,'Min/Max rhodrysig =',minval(rhodrysig),maxval(rhodrysig)
      print *,'Min/Max zsig =',minval(hzsig),maxval(hzsig)
   endif
   deallocate(pcpmr_sig,rhodrysig)
 endif ! large_ngrid

! Generate derived reflectivity fields.
 if(.true.)then ! always run this part
   if (verbose) then
      print*,' '
      print *,'Calling lfm_reflectivity.'
      print *,'Min/Max hrainmr_sig =',minval(hrainmr_sig),maxval(hrainmr_sig)
      print *,'Min/Max hcldicemr_sig =',minval(hcldicemr_sig),maxval(hcldicemr_sig)
      print *,'Min/Max hsnowmr_sig =',minval(hsnowmr_sig),maxval(hsnowmr_sig)
      print *,'Min/Max hgraupelmr_sig =',minval(hgraupelmr_sig),maxval(hgraupelmr_sig)
   endif
   if (minval(hrainmr_sig) < rmsg .and. minval(hcldicemr_sig) < rmsg .and.  &
       minval(hsnowmr_sig) < rmsg .and. minval(hgraupelmr_sig) < rmsg) then
      allocate(tvsig(lx,ly,nz),rhomoistsig(lx,ly,nz))
      tvsig=htsig*(1.+0.61*hmrsig)
      rhomoistsig=hpsig/(r*tvsig)
      call lfm_reflectivity(lx,ly,nz,rhomoistsig,hzsig  &
                           ,hrainmr_sig,hcldicemr_sig,hsnowmr_sig,hgraupelmr_sig  &
                           ,hrefl_sig,hzdr_sig,hldr_sig,htsig,hpsig,hmrsig        &
                           ,max_refl,echo_tops)

!     Determine low level reflectivity
      if(.false.)then ! use lowest sigma level (sfc)
          refl_sfc(:,:)=hrefl_sig(:,:,1)
      elseif(large_ngrid)then ! assume a certain sigma level is 1km AGL
          refl_sfc(:,:)=hrefl_sig(:,:,8)
      else ! interpolate to 1km AGL
          do j=1,ly
          do i=1,lx
              ht_1km = zsfc(i,j) + 1000.
              do k = 1,nz-1
                  if(hzsig(i,j,k) .le. ht_1km .AND. hzsig(i,j,k+1) .ge. ht_1km)then
                      frack = (ht_1km - hzsig(i,j,k)) / (hzsig(i,j,k+1) - hzsig(i,j,k))
                      klow = k
                      khigh = klow+1
                      refl_sfc(i,j)=hrefl_sig(i,j,klow) * (1. - frack) + hrefl_sig(i,j,khigh) * frack
                  endif
              enddo ! k
          enddo ! i
          enddo ! j
      endif

!     Run Advection routine on reflectivity
      if(fcsttime .le. i4_adv_pcp)then
          I4_elapsed = ishow_timer()

          ext = 'lmr'
          var_2d = 'R'
          call get_laps_2d(laps_reftime,ext,var_2d,units_2d,comment_2d,lx,ly,max_refl,istatus)
          if(istatus .eq. 1)then
            write(6,*)' Obtained initial time max_refl field from lmr analysis LMR file'
          else
            write(6,*)' Warning: could not find initial time max_refl field from lmr analysis LMR file'
          endif      

          call mean_wind_bunkers(husig,hvsig,hzsig(:,:,1),lx,ly,nz &     ! I
                          ,hzsig                                   &     ! I
                          ,umean,vmean,ustorm,vstorm,status)             ! O
          write(6,*)' Advect Surface Reflectivity ',fcsttime,i4_adv_pcp
          call advect(umean,vmean,refl_sfc,array_buf,ngrid_spacingx & 
                     ,lx,ly,array_out,float(fcsttime),1.0,lon,rmsg)
          refl_sfc = array_out
          write(6,*)' Advect Max Reflectivity ',fcsttime,i4_adv_pcp
          call advect(umean,vmean,max_refl,array_buf,ngrid_spacingx & 
                     ,lx,ly,array_out,float(fcsttime),1.0,lon,rmsg)
          max_refl = array_out
      else
          write(6,*)' Using model output reflectivity (no advection)',fcsttime,i4_adv_pcp
      endif

      if(fcsttime .le. i4_adv_cld)then
          I4_elapsed = ishow_timer()
          write(6,*)' Advect 3D Clouds ',fcsttime,i4_adv_cld
          do k = 1,nz
              call advect(husig(:,:,k),hvsig(:,:,k),hcldliqmr_sig(:,:,k),array_buf,ngrid_spacingx & 
                         ,lx,ly,array_out,float(fcsttime),1.0,lon,rmsg)
              hcldliqmr_sig(:,:,k) = array_out(:,:)

              call advect(husig(:,:,k),hvsig(:,:,k),hcldicemr_sig(:,:,k),array_buf,ngrid_spacingx & 
                         ,lx,ly,array_out,float(fcsttime),1.0,lon,rmsg)
              hcldicemr_sig(:,:,k) = array_out(:,:)
          enddo ! k
      else
          write(6,*)' Using model output clouds (no advection)',fcsttime,i4_adv_cld
      endif

      I4_elapsed = ishow_timer()

!     Conversion of microphysical cloud mixing ratios to concentrations
      hcldliqmr_sig(:,:,:) = hcldliqmr_sig(:,:,:) * rhomoistsig(:,:,:)
      hcldicemr_sig(:,:,:) = hcldicemr_sig(:,:,:) * rhomoistsig(:,:,:)

      deallocate(tvsig,rhomoistsig)
      if (verbose) then
        print*,' '
        print *,'After lfm_reflectivity.'
        print *,'Min/Max hrefl_sig =',minval(hrefl_sig),maxval(hrefl_sig)
        print *,'Min/Max hzdr_sig =',minval(hzdr_sig),maxval(hzdr_sig)
        print *,'Min/Max hldr_sig =',minval(hldr_sig),maxval(hldr_sig)
        print *,'Min/Max refl_sfc =',minval(refl_sfc),maxval(refl_sfc)
      endif
   else
      hrefl_sig=rmsg
      hzdr_sig=rmsg
      hldr_sig=rmsg
      max_refl=rmsg
      echo_tops=rmsg
   endif


! Do some QC on cloud species fields.

   where(cldliqmr_prs < zero_thresh) cldliqmr_prs=0.
   where(cldicemr_prs < zero_thresh) cldicemr_prs=0.
   where(snowmr_prs < zero_thresh) snowmr_prs=0.
   where(rainmr_prs < zero_thresh) rainmr_prs=0.
   where(graupelmr_prs < zero_thresh) graupelmr_prs=0.

! Generate surface and ua precip type.

   if (verbose) then
      print*,' '
      print *,'Calling lfm_ua_pcptype.'
      print *,'Min/Max hrainmr_sig =',minval(hrainmr_sig),maxval(hrainmr_sig)
      print *,'Min/Max hsnowmr_sig =',minval(hsnowmr_sig),maxval(hsnowmr_sig)
      print *,'Min/Max hgraupelmr_sig =',minval(hgraupelmr_sig),maxval(hgraupelmr_sig)
   endif
   if (minval(hrainmr_sig) < rmsg .and.                            &
       minval(hsnowmr_sig) < rmsg .and. minval(hgraupelmr_sig) < rmsg) then

      call lfm_sfc_pcptype(lx,ly,nz,hrainmr_sig,hsnowmr_sig,hgraupelmr_sig  &
                          ,pcptype_sfc)
      if (verbose) then
        print*,' '
        print *,'After lfm_sfc_pcptype.'
        print *,'Min/Max pcptype_sfc =',minval(pcptype_sfc),maxval(pcptype_sfc)
      endif

      if(large_pgrid .eqv. .false.)then
        call lfm_ua_pcptype(lx,ly,lz,rainmr_prs,snowmr_prs,graupelmr_prs  &
                           ,pcptype_prs)
        if (verbose) then
          print*,' '
          print *,'After lfm_ua_pcptype.'
          print *,'Min/Max pcptype_prs =',minval(pcptype_prs),maxval(pcptype_prs)
        endif
      endif
 
   else
      pcptype_sfc=rmsg
      pcptype_prs=rmsg
   endif

!  Convert mr to concentration for microphysical species?
!  To do this we'd need rho on the pressure grid. Otherwise we can convert
!  on the model sig grid using rhosig, and this is now done at the end of 
!  'lfm_reflectivity' (and right after returning)

 endif ! true

endif ! make_micro

I4_elapsed = ishow_timer()

! Precip fields.

do k=1,lz
   if (lprs(k) == 100000) k1000=k
   if (lprs(k) == 85000) k850=k
   if (lprs(k) == 70000) k700=k
   if (lprs(k) == 50000) k500=k
enddo

allocate(pcp_init(lx,ly))
allocate(pcp_06(lx,ly))
call fill_precip(pcp_init,pcp_06,snow_tot)

! fcsttime is in seconds
! change nmm, nmm zeroes accumulated precipitation every interval period

if (fcsttime > 0) then
   if (trim(mtype) == 'nmm') then
       pcp_inc=pcp_tot
       print*, 'pcp_inc from pcp_tot: ',trim(mtype)
   else
     pcp_inc=pcp_tot-pcp_init
     print*, 'pcp_inc from pcp_init: ',trim(mtype)
     print*, 'fcsttime: ',fcsttime
     print*,maxval(pcp_inc),maxval(pcp_tot),maxval(pcp_init)
   endif
   where(pcp_inc < .0001) pcp_inc=0.
   where(pcp_tot < .0001) pcp_tot=0.

   if(.not. large_ngrid)then
      print*,'Min/Max snow tot (previous) = ',minval(snow_tot),maxval(snow_tot)

!     allocate(fallen_precip_type(lx,ly))
!     call wintprec(htsig,hzsig,zprs,psfc,tsfc,tdsfc,zsfc,pcp_inc,lx,ly  &
!               ,nz,lz,k700,k850,k1000,pcp_inc,fallen_precip_type)
!     call snowfall(htsig,tsfc,pcp_inc,fallen_precip_type,lx,ly,nz,snow_inc,snow_tot)
      call snowfall(htsig,tsfc,pcp_inc,pcptype_sfc,lx,ly,nz,snow_inc,snow_tot)
!     deallocate(fallen_precip_type)
   endif

else ! no incremental precip at time zero
   pcp_inc = r_missing_data

endif

if (trim(mtype) == 'nmm') then
 print*, 'pcp_inc update: ',trim(mtype)
 pcp_06 = pcp_tot
 pcp_tot = 0.
 pcp_tot = pcp_init + pcp_inc
else
 print*, 'skip pcp_inc update: ',trim(mtype)
!pcp_tot = pcp_tot + pcp_inc
endif

! Write intermediate precip file for future use.

open(1,file=trim(filename0),status='unknown',form='unformatted')
if (trim(mtype) == 'nmm') then
  write(1) pcp_tot,pcp_06,snow_tot
  print*,'write intermediate pcp file for model: ',trim(mtype)
else
  write(1) pcp_tot,snow_tot
endif
close(1)
if (verbose) then
   print*,' '
   print*,'Writing new precip data to: ',trim(filename0)
endif
deallocate(pcp_init)
deallocate(pcp_06)

! Adjust precip/snow fields with a constant (for now) bias correction
const_pcpadj = 1.00 ! 0.65
where(pcp_inc .ne. r_missing_data)pcp_inc = pcp_inc * const_pcpadj
where(pcp_tot .ne. r_missing_data)pcp_tot = pcp_tot * const_pcpadj
where(snow_inc .ne. r_missing_data)snow_inc = snow_inc * const_pcpadj
where(snow_tot .ne. r_missing_data)snow_tot = snow_tot * const_pcpadj

! Temporary variables needed to derive some fields.

allocate(hrhsig(lx,ly,nz),hthetasig(lx,ly,nz),hthetaesig(lx,ly,nz),hzsigf(lx,ly,nz+1))
do k=1,nz
do j=1,ly
do i=1,lx
   hrhsig(i,j,k)=min(1.,relhum(htsig(i,j,k),hmrsig(i,j,k),hpsig(i,j,k)))
   hthetasig(i,j,k)=potential_temp(htsig(i,j,k),hpsig(i,j,k))
   hthetaesig(i,j,k)=eq_potential_temp(htsig(i,j,k),hpsig(i,j,k)  &
                                      ,hmrsig(i,j,k),hrhsig(i,j,k))
enddo
enddo
enddo
do k=2,nz
   hzsigf(:,:,k)=(hzsig(:,:,k-1)+hzsig(:,:,k))*0.5
enddo
hzsigf(:,:,1)=zsfc
hzsigf(:,:,nz+1)=2.*hzsig(:,:,nz)-hzsigf(:,:,nz)

! Generate pbl height if not available from model.

if (maxval(pblhgt) > 999999.) then
   if (verbose) then
      print*,' '  
      print*,'Generating PBL height using theta and surface temp.'
   endif

   do j=1,ly
   do i=1,lx
      thetasfc(i,j)=potential_temp(tsfc(i,j),psfc(i,j))
   enddo
   enddo
   call model_pblhgt(hthetasig,thetasfc,hpsig,hzsig,zsfc,lx,ly,nz,pblhgt)
endif
where (pblhgt < 0.) pblhgt=0.

I4_elapsed = ishow_timer()

if(.not. large_ngrid)then

! Helicity, cape, cin, LI.

  write(6,*)' Calculating helicity, stability indices'
  call helicity(husig,hvsig,hzsig,usfc,vsfc,zsfc,lx,ly,nz,srhel)

  I4_elapsed = ishow_timer()

  call updraft_helicity(husig,hvsig,hwsig,hzsig,hzsigf,zsfc,llat,llon,lx,ly,nz,uhel)
  print *,'Min/Max uhel =',minval(uhel),maxval(uhel)

  call bulk_shear(husig,hvsig,hzsig,hzsigf,zsfc,llat,llon,lx,ly,nz,bshr)

  I4_elapsed = ishow_timer()

  call capecin(hpsig*0.01,htsig,hthetaesig,hthetasig,hrhsig  &
              ,hzsigf,tprs,liftedind,cape,cin,k500,lx,ly,nz,lz)
 
! supercell composite parameter
  do i = 1,lx
  do j = 1,ly
    if(bshr(i,j) .ge. 20.)then
      arg = 1.0
    elseif(bshr(i,j) .ge. 10.)then
      arg = bshr(i,j)/20.
    else
      arg = 0.0
    endif
    scp(i,j) = cape(i,j)/1000. * (srhel(i,j)/50.) * arg
  enddo ! j
  enddo ! i
endif ! large_ngrid

deallocate(hthetasig,hthetaesig,hzsigf)

I4_elapsed = ishow_timer()

! Height of wet-bulb = 0C.

write(6,*)' Calculating wet bulb zero height'
allocate(htdsig(lx,ly,nz))
if(.not. large_ngrid)then
  htdsig=htsig/((-rvolv*alog(hrhsig)*htsig)+1.0)

  call height_tw(hpsig,hzsig,htsig,htdsig,psfc,zsfc,tsfc,tdsfc,mrsfc,pmsl  &
                ,0.0,ztw0,lx,ly,nz)

  I4_elapsed = ishow_timer()

  call height_tw(hpsig,hzsig,htsig,htdsig,psfc,zsfc,tsfc,tdsfc,mrsfc,pmsl  &
                ,1.3,ztw1,lx,ly,nz)

  I4_elapsed = ishow_timer()

  aglhgt = 80.
  call wind_agl(husig,hvsig,hzsig,aglhgt,zsfc,lx,ly,nz  &
               ,u80,v80,r_missing_data)

endif ! large_ngrid

! Visibility.

if(.true.)then
! a1 = 75.; b1 = 37.5; c1 = 1.5; d1 = 5.357; e1 = 1.2

  clwc2alpha = 1.5 / (rholiq  * reff_clwc)
  cice2alpha = 1.5 / (rholiq  * reff_cice)
  rain2alpha = 1.5 / (rholiq  * reff_rain)
  snow2alpha = 1.5 / (rhosnow * reff_snow)
  pice2alpha = 1.5 / (rhograupel * reff_graupel)

  a1 = clwc2alpha
  b1 = cice2alpha
  c1 = rain2alpha
  d1 = snow2alpha
  e1 = pice2alpha

! Note that hydrometeor mixing ratios have previously been converted to content
  alpha(:,:) = a1 * hcldliqmr_sig(:,:,1) + b1 * hcldicemr_sig(:,:,1) &
             + c1 * hrainmr_sig(:,:,1) + d1 * hsnowmr_sig(:,:,1) + e1 * hgraupelmr_sig(:,:,1)
! meteorological visibility, tau ~= 2.8, with lower bound added on alpha
  visibility(:,:) = 2.8 / max(alpha(:,:),.00004) 
else
  visibility=6000000.*(tsfc-tdsfc)/(rhsfc**1.75)  ! in meters
  where(visibility > 99990.) visibility = 99990.
endif

! Compute heat index if temp is above 80F (300K).

do j=1,ly
do i=1,lx
   if (tsfc(i,j) >= 300.) then
      heatind(i,j)=heatindex(tsfc(i,j),rhsfc(i,j))
   else
      heatind(i,j)=tsfc(i,j)
   endif
enddo
enddo

! Fire weather indices.

if (make_firewx .and. .not. large_ngrid) then
   call ventilation(husig,hvsig,hzsig,pblhgt,zsfc,lx,ly,nz,upbl,vpbl,vnt_index)

! Mid-level Haines Index.

   call haines_layer(hpsig,htsig,htdsig,ham_index,lx,ly,nz,850.,700.)

! High-level Haines Index.

   call haines_layer(hpsig,htsig,htdsig,hah_index,lx,ly,nz,700.,500.)

! Fosberg FWI.

   call fosberg_fwi(tsfc,rhsfc,usfc,vsfc,lx,ly,fwi_index)

! Kelsch FWX.
   do i = 1,lx
   do j = 1,ly
       windspeed(i,j)=sqrt(usfc(i,j)**2+vsfc(i,j)**2)
   enddo ! j
   enddo ! i
   ismoist = 0
   isnow = 0
   call lp_fire_danger (lx,ly,rhsfc,tsfc,windspeed, &
                        soil_moist,snow_cover,zsfc,ldf,  &
                        ismoist,isnow,fwx_index,istatus)

endif

if(allocated (hrhsig)) deallocate(hrhsig)
if(allocated (htdsig)) deallocate(htdsig)

!cj Compute Omega, added 6/21/2007
if(.not. large_pgrid)then
  do k=1,lz
  do j=1,ly
  do i=1,lx
   tvprs = tprs(i,j,k) * (1. + 0.61*shprs(i,j,k))
!  YHL  omprs(i,j,k) = -(lprsl(k)*wprs(i,j,k)*grav)/(r*tvprs) * 100. ! Pa/Mb
   omprs(i,j,k) = -(lprs(k)*wprs(i,j,k)*grav)/(r*tvprs)  ! changed by Steve. A, Huiling Yuan 9/11/2009
  enddo
  enddo
  enddo
endif

if(fcsttime .le. i4_adv_cld)then
  call mean_wind_bunkers(husig,hvsig,hzsig(:,:,1),lx,ly,nz &     ! I
                        ,hzsig                                   &     ! I
                        ,umean,vmean,ustorm,vstorm,status)             ! O
  ext = 'lcv'
  var_2d = 'SWI'
  call get_laps_2d(laps_reftime,ext,var_2d,units_2d,comment_2d,lx,ly,swdown,istatus)
  if(istatus .eq. 1)then
    write(6,*)' Obtained initial time swdown field from swi analysis LCV file'
  else
    write(6,*)' Warning: could not find initial time swdown field from swi analysis LCV file'
  endif      
  write(6,*)' Advect Solar Radiation ',fcsttime,i4_adv_cld
  call advect(umean,vmean,swdown,array_buf,ngrid_spacingx & 
             ,lx,ly,array_out,float(fcsttime),1.0,lon,rmsg)
  write(6,*)' Range of advection input is ',minval(swdown),maxval(swdown)
  write(6,*)' Range of initial advection output is ',minval(array_out),maxval(array_out)

  write(6,*)' Calling get_ghi_ratio for times ',laps_reftime,laps_valtime
  call get_ghi_ratio(laps_reftime,laps_valtime,lat,lon,lx,ly,ghi_ratio)
  write(6,*)' Range of ghi_ratio is ',minval(ghi_ratio),maxval(ghi_ratio)

  swdown = array_out * ghi_ratio
  write(6,*)' Range of final advection output is ',minval(swdown),maxval(swdown)
else
  if (fcsttime .eq. 0) then ! obtain initial time swdown field from swi analysis LCV files
    ext = 'lcv'
    var_2d = 'SWI'
    call get_laps_2d(laps_reftime,ext,var_2d,units_2d,comment_2d,lx,ly,swdown,istatus)
    if(istatus .eq. 1)then
      write(6,*)' Obtained initial time swdown field from swi analysis LCV file'
    else
      write(6,*)' Warning: could not find initial time swdown field from swi analysis LCV file'
    endif    
  else
    if(.true.)then ! post-process swdown field  
       do j=1,ly
       do i=1,lx
          coeff = 0.93 + (cldamt(i,j) * 0.5)
          swdown(i,j) = swdown(i,j) * coeff
       enddo ! i
       enddo ! j      
    endif
  endif
endif  

if (verbose .and. .not. large_ngrid) then
   print*,' '
   print*,'Min/Max tsfc      = ',minval(tsfc),maxval(tsfc)
   print*,'Min/Max rhsfc     = ',minval(rhsfc),maxval(rhsfc)
   print*,'Min/Max tdsfc     = ',minval(tdsfc),maxval(tdsfc)
   print*,'Min/Max thetasfc  = ',minval(thetasfc),maxval(thetasfc)
   print*,'Min/Max thetaesfc = ',minval(thetaesfc),maxval(thetaesfc)
   print*,'Min/Max redp      = ',minval(redp),maxval(redp)
   print*,'Min/Max pmsl      = ',minval(pmsl),maxval(pmsl)
   print*,'Min/Max ztw0      = ',minval(ztw0),maxval(ztw0)
   print*,'Min/Max ztw1      = ',minval(ztw1),maxval(ztw1)
   print*,'Min/Max pblhgt    = ',minval(pblhgt),maxval(pblhgt)
   print*,'Min/Max cldbase   = ',minval(cldbase),maxval(cldbase)
   print*,'Min/Max cldtop    = ',minval(cldtop),maxval(cldtop)
   print*,'Min/Max ceiling   = ',minval(ceiling),maxval(ceiling)
   print*,'Min/Max cldamt    = ',minval(cldamt),maxval(cldamt)
   print*,'Min/Max intliqwat = ',minval(intliqwater),maxval(intliqwater)
   print*,'Min/Max intcldice = ',minval(intcldice),maxval(intcldice)
   print*,'Min/Max totpcpwat = ',minval(totpcpwater),maxval(totpcpwater)
   print*,'Min/Max intrain   = ',minval(intrain),maxval(intrain)
   print*,'Min/Max intsnow   = ',minval(intsnow),maxval(intsnow)
   print*,'Min/Max intgraupel= ',minval(intgraupel),maxval(intgraupel)
   print*,'Min/Max simvis    = ',minval(simvis),maxval(simvis)
   print*,'Min/Max max refl  = ',minval(max_refl),maxval(max_refl)
   print*,'Min/Max echo tops = ',minval(echo_tops),maxval(echo_tops)
   print*,'Min/Max refl sfc  = ',minval(refl_sfc),maxval(refl_sfc)
   print*,'Min/Max pcptypsfc = ',minval(pcptype_sfc),maxval(pcptype_sfc)
   print*,'Min/Max pcp inc   = ',minval(pcp_inc),maxval(pcp_inc)
   print*,'Min/Max snow inc  = ',minval(snow_inc),maxval(snow_inc)
   print*,'Min/Max pcp tot   = ',minval(pcp_tot),maxval(pcp_tot)
   print*,'Min/Max snow tot  = ',minval(snow_tot),maxval(snow_tot)
   print*,'Min/Max helicity  = ',minval(srhel),maxval(srhel)
   print*,'Min/Max cape      = ',minval(cape),maxval(cape)
   print*,'Min/Max cin       = ',minval(cin),maxval(cin)
   print*,'Min/Max li        = ',minval(liftedind),maxval(liftedind)
   print*,'Min/Max vis       = ',minval(visibility),maxval(visibility)
   print*,'Min/Max heat ind  = ',minval(heatind),maxval(heatind)
   print*,'Min/Max lw out    = ',minval(lwout),maxval(lwout)
   print*,'Min/Max sw out    = ',minval(swout),maxval(swout)
   print*,'Min/Max lw down   = ',minval(lwdown),maxval(lwdown)
   print*,'Min/Max sw down   = ',minval(swdown),maxval(swdown)
   print*,'Min/Max sh flux   = ',minval(shflux),maxval(shflux)
   print*,'Min/Max lh flux   = ',minval(lhflux),maxval(lhflux)
   if (make_firewx) then
      print*,'Min/Max upbl      = ',minval(upbl),maxval(upbl)
      print*,'Min/Max vpbl      = ',minval(vpbl),maxval(vpbl)
      print*,'Min/Max vent indx = ',minval(vnt_index),maxval(vnt_index)
      print*,'Min/Max md haines = ',minval(ham_index),maxval(ham_index)
      print*,'Min/Max up haines = ',minval(hah_index),maxval(hah_index)
      print*,'Min/Max fosberg   = ',minval(fwi_index),maxval(fwi_index)
      print*,'Min/Max laps/kelsch = ',minval(fwx_index),maxval(fwx_index)
   endif
endif

I4_elapsed = ishow_timer()

!Beka moisture flux

        call get_grid_spacing_array(lat,lon,lx,ly,dx,dy)

!       print *,'Min/Max zsig =',minval(hzsig),maxval(hzsig)
!       write(6,*)'lx,ly,lz',lx,ly,lz

!       do i = 1,lx
!       do j = 1,ly
!           if(hzsig(i,j,1) .gt. hzsig(i,j,lz))then
!               write(6,*)'hzsig QC issue ',i,j 
!           endif
!       enddo ! j
!       enddo ! i

if(.not. large_ngrid)then

        call up_mflux(lx,ly,nz,zsfc,ldf,dx,dy                 &
                     ,husig,hvsig,totpcpwater,upflux           &
                     ,hzsig,rmsg)              

!       write(*,*)'upflux/totpcpwater fields'

!       write(*,*)upflux,totpcpwater

	write(*,*)'beka beka beka'

endif

if (fcsttime .eq. 0) then ! obtain initial time bt11u field from s8a analysis LCV files
  ext = 'lcv'
  var_2d = 'S8A'
  call get_laps_2d(laps_reftime,ext,var_2d,units_2d,comment_2d,lx,ly,bt11u,istatus)
  if(istatus .eq. 1)then
    write(6,*)' Obtained initial time bt11u field from s8a analysis LCV file'
  else
    write(6,*)' Warning: could not find initial time swdown field from s8a analysis LCV file'
  endif      

else ! Calculate brightness temperature
  stefan_boltzmann = 5.67e-8
  if(.false.)then ! First method                                 
    eff_emissivity = 0.6
    bt11u(:,:) = ( (lwout(:,:)/eff_emissivity) / stefan_boltzmann) ** 0.25 
  else          ! Second method: adapted from Ohring, George, Arnold Gruber, Robert Ellingson, 1984: 
                ! Satellite Determinations of the Relationship between Total Longwave Radiation Flux 
                ! and Infrared Window Radiance. J. Climate Appl. Meteor., 23, 416-425.
                ! Based on a look at graphs on this and another paper a quadratic fit is being used    
                ! using these three points: 200,200 240,255 and 280,320
    do j=1,ly
    do i=1,lx
        bt_flux_equiv = (lwout(i,j) / stefan_boltzmann) ** 0.25 
        bt11u(i,j) = 255. + 1.5 * (bt_flux_equiv - 240.) + .003125 * (bt_flux_equiv - 240.)**2
    enddo ! i
    enddo ! j
  endif

endif

return
end

!===============================================================================

subroutine fill_precip(pcp_init,pcp_06,snow_init)

use lfmgrid

implicit none

integer :: nc, itry
real, dimension(lx,ly) :: pcp_init,pcp_06,snow_init
character(len=10) :: atime
logical :: back,there

! Read previous accumulated precip data from an intermediate file that
!  is created in the model run directory.

back=.true.
nc=index(filename,'/',back)

if (fcsttime > 0) then
   if (nc > 0) filename0=filename(1:nc)//'/'//domnum_fstr
   write(atime,'(i6.6)') max(0,fcsttime-precip_dt)
   filename0=trim(filename0)//'_'//trim(adjustl(atime))//'.pcp'
   itry = 0

80 inquire(file=trim(filename0),exist=there)
   if (there) then
90    if (verbose) then
         print*,' '
         print*,'Reading previous precip data from: ',trim(filename0)
      endif
      open(1,file=trim(filename0),status='old',form='unformatted')
      if (trim(mtype) == 'nmm') then
        read(1,end=101,err=101) pcp_init,pcp_06,snow_init
        print *,'reading intermediate file for model: ',trim(mtype)
      else
        read(1,end=101,err=101) pcp_init,snow_init
        pcp_06=0.
      endif
      goto 102

101   print*,'  ERROR reading previous precip'

      itry = itry + 1
      if(itry .le. 5)then
          print*,'  trying once again after waiting'          
          call sleep(60)
          goto 90
      endif

      print*,'  Prior precip set to zero.'
      pcp_init=0.
      snow_init=0.

102   close(1)

   else ! not there
      print*,'Could not find previous precip file: ',trim(filename0)

      itry = itry + 1
      if(itry .le. 5)then
          print*,'  trying once again after waiting'          
          call sleep(60)
          goto 80
      endif

      print*,'  Prior precip set to zero.'
      pcp_init=0.
      snow_init=0.
   endif

else ! fcsttime = 0
   pcp_init=0.
   snow_init=0.
   pcp_inc=0.
   snow_inc=0.
   pcp_tot=0.
   snow_tot=0.

endif

! Create new intermediate filename for current accumulated precip data.

if (nc > 0) filename0=filename(1:nc)//'/'//domnum_fstr
write(atime,'(i6.6)') fcsttime
filename0=trim(filename0)//'_'//trim(adjustl(atime))//'.pcp'

return
end

!===============================================================================

subroutine model_pblhgt(theta,thsfc,psig,zsig,topo,nx,ny,nz,pblhgt)

!  Subroutine to estimate height AGL in meters of PBL from native
!  coordinate model data.  Adapted from the LAPS routine for 
!  terrain-following model coordinates.

implicit none

integer :: nx,ny,nz,i,j,k,ktop
real :: thresh_k,topwgt,botwgt
real, dimension(nx,ny,nz) :: theta(nx,ny,nz),psig,zsig
real, dimension(nx,ny) :: thsfc,topo,pblhgt
logical :: found_pbl_top

do j=1,ny
do i=1,nx

! Compute threshold value that theta needs to exceed
! to be above PBL.  We use surface theta plus an
! additional 3K for slop.

   thresh_k = thsfc(i,j) + 3.0  

! Now begin at the bottom and work our way up until we find 
! the first level with a theta exceeding the threshold.

   found_pbl_top = .false.
   do k=1,nz
 
      if (theta(i,j,k) >= thresh_k) then
         ktop = k
         found_pbl_top = .true.
         exit 
      endif
   enddo 

! If we did not find a good PBL, set PBL to first level
! and print out some diagnostics.

   if (.not. found_pbl_top) then
      print*,'PBL height not found at i/j = ',i,j
      print*,'Surface Theta = ',thsfc(i,j)
      print*,'Theta in the column:'
      print*,'Pressure Height  Theta'
      print*,'-------- ------- --------'
      do k=1,nz
         print '(F8.0,F8.0,F8.2)',psig(i,j,k),zsig(i,j,k),theta(i,j,k)
      enddo
      ktop = 1
      pblhgt(i,j)=zsig(i,j,1)-topo(i,j)
   else 
  
! We found the top k-level bounding the PBL so interpolate
! to the actual level.
 
      if (ktop == 1) then
         pblhgt(i,j)=zsig(i,j,1)-topo(i,j)
      else
! Interpolate to get height at thresh_k.
         botwgt=((theta(i,j,ktop)-thresh_k)  &
                /(theta(i,j,ktop)-theta(i,j,ktop-1)))
         topwgt=1.0 - botwgt
         pblhgt(i,j)=botwgt*zsig(i,j,ktop-1)  &
                    +topwgt*zsig(i,j,ktop)-topo(i,j)
      endif
   endif
enddo
enddo

return
end

!===============================================================================

subroutine lfmclouds(nx,ny,nz,grid_spacing,cldliqmr,cldicemr,snowmr  &
                    ,height,topo,cldbase,cldtop,ceiling,cldamt)

! This routine is used to compute cloud ceiling (AGL), cloud base
! height (MSL), cloud top height (MSL), and coverage (fraction) given
! mixing ratios of the various microphysical species.  

! Adapted from the USAF Weather Agency MM5 post processer.
!  Brent Shaw, NOAA Forecast Systems Lab, Dec 2000

implicit none

integer :: nx,ny,nz,i,j,k
real :: grid_spacing,icethresh,liqthresh,snowthresh
real, dimension(nx,ny) :: topo,cldbase,cldtop,ceiling,cldamt
real, dimension(nx,ny,nz) :: cldliqmr,cldicemr,snowmr,height    

! Set thresholds for cloud determination based on grid
! resolution.  Kind of hokey, but will get us by for now.

if (grid_spacing <= 10000) then
   icethresh=0.000005
   snowthresh=0.000003
   liqthresh=0.000003
else
   icethresh=0.000005
   snowthresh=0.000025
   liqthresh=0.000025
endif

! Loop through using these thresholds to determine cloudiness.

do j=1,ny
do i=1,nx
   do k=1,nz
      if ((cldliqmr(i,j,k) >= liqthresh) .or.  &
          (cldicemr(i,j,k) >= icethresh) .or.  &
          (snowmr(i,j,k)   >= snowthresh) ) then
         cldbase(i,j)=height(i,j,k)
         ceiling(i,j)=height(i,j,k)-topo(i,j)
              
! For now all we can do is use cldamt as a yes/no.
! We should look at coming up with a fraction function.

         cldamt(i,j) = 1.0
         exit
      endif
   enddo
   do k=nz,1,-1
      if ((cldliqmr(i,j,k) >= liqthresh) .or.  &
          (cldicemr(i,j,k) >= icethresh) .or.  &
          (snowmr(i,j,k)   >= snowthresh)) then
         cldtop(i,j)=height(i,j,k)
         exit
      endif
   enddo
enddo
enddo    
    
return
end

!===============================================================================

subroutine lfm_integrated_liquid(nx,ny,nz,cond_mr,cice_mr,vapor_mr,rain_mr,snow_mr,graupel_mr,rho,height,topo  &
                                ,intliqwater,intcldice,totpcpwater,intrain,intsnow,intgraupel,rmsg)

! Computes integrated liquid water and total precip. water in a column.  
!  Adapted from USAF Weather Agency MM5 Post Processor
!  Brent Shaw, NOAA Forecast Systems Lab, Dec 2000

implicit none
  
integer :: nx,ny,nz,i,j,k 
real :: height_top,height_bot,dz,rmsg
real, dimension(nx,ny) :: newtopo,intliqwater,totpcpwater,intcldice  ! M
real, dimension(nx,ny) :: topo,intrain,intsnow,intgraupel         ! M
real, dimension(nx,ny,nz) :: cond_mr,vapor_mr,cice_mr             ! KG/KG
real, dimension(nx,ny,nz) :: rain_mr,snow_mr,graupel_mr           ! KG/KG
real, dimension(nx,ny,nz) :: rho                                  ! KG/M**3
real, dimension(nx,ny,nz) :: height                               ! M
real rho_h2o                                                      ! KG/M**3

rho_h2o = 1000.                                                   ! KG/M**3

do j=1,ny
do i=1,nx
   intliqwater(i,j)=0.0
   totpcpwater(i,j)=0.0
   intcldice(i,j)=0.0
   intrain(i,j)=0.0
   intsnow(i,j)=0.0
   intgraupel(i,j)=0.0
   do k=1,nz
! Compute layer thickness
      if (k == 1) then
         height_bot=topo(i,j)
         height_top=0.5*(height(i,j,1)+height(i,j,2))
      elseif (k == nz) then
         height_bot=0.5*(height(i,j,nz-1)+height(i,j,nz))
         height_top=2*height(i,j,nz)-height_bot
      else
         height_bot=0.5*(height(i,j,k-1)+height(i,j,k))
         height_top=0.5*(height(i,j,k)+height(i,j,k+1))
      endif
      dz=height_top-height_bot
      if(cond_mr(i,j,k) < rmsg)intliqwater(i,j)=intliqwater(i,j)+cond_mr(i,j,k) *(rho(i,j,k)/rho_h2o)*dz  ! meters
      if(cice_mr(i,j,k) < rmsg)intcldice(i,j)  =intcldice(i,j)  +cice_mr(i,j,k) *(rho(i,j,k)/rho_h2o)*dz  ! meters
      if(rain_mr(i,j,k) < rmsg)intrain(i,j)    =intrain(i,j)    +rain_mr(i,j,k) *(rho(i,j,k)/rho_h2o)*dz  ! meters
      if(snow_mr(i,j,k) < rmsg)intsnow(i,j)    =intsnow(i,j)    +snow_mr(i,j,k) *(rho(i,j,k)/rho_h2o)*dz  ! meters
      if(graupel_mr(i,j,k) < rmsg)intgraupel(i,j) =intgraupel(i,j)    +graupel_mr(i,j,k) *(rho(i,j,k)/rho_h2o)*dz  ! meters
      totpcpwater(i,j)=totpcpwater(i,j)                         +vapor_mr(i,j,k)*(rho(i,j,k)/rho_h2o)*dz  ! meters
   enddo
enddo
enddo

return
end 

!===============================================================================

subroutine lfm_reflectivity(lx,ly,nz,rho,hgt,rainmr,icemr,snowmr,graupelmr  &
                           ,refl,zdr,ldr,tmp,p,qv,max_refl,echo_tops)

use lfmgrid, ONLY: c_m2z, large_ngrid
use src_versuch, ONLY: versuch

! Subroutine to compute estimated radar reflectivity from
! the precipitation mixing ratios.  Will also return 
! column max reflectivity and echo tops.  The estimation
! is done using formulae from Kessler (1969) and 
! Rogers and Yau (1989).  

! Adapted from USAF Weather Agency routine.  
!  Brent Shaw, NOAA Forecast System Lab, Dec 2000
  
implicit none

integer :: lx,ly,nz,i,j,k
integer :: in0r,in0s,in0g,iliqskin
real, parameter :: svnfrth=7.0/4.0,max_top_thresh=5.0
real :: w
real, dimension(lx,ly) :: max_refl,echo_tops
real, dimension(lx,ly,nz) :: rho,rainmr,icemr,snowmr,graupelmr,refl
real, dimension(lx,ly,*)  :: hgt
real, dimension(lx,ly,nz) :: zdr,ldr
real, dimension(lx,ly,nz) :: tmp,p,qv

! Local arrays

real, dimension(lx,ly,nz) :: zhhx,ahhx,zhvx,zvhx,zvvx    ! versuch outputs
real, dimension(lx,ly,nz) :: delahvx,kkdp,rhohvx         ! versuch outputs

real :: zvvxmax,zhhxmax

max_refl=0.0
echo_tops=1.0e37
if (c_m2z .ne. 'wrf') then
    refl=0.0
else
    refl = max(refl,-10.)
    where(refl .lt. 0. .and. refl .gt. -10.)refl = 0.
endif
zdr=0.0
ldr=0.0

print *,'c_m2z = ',c_m2z   

if(c_m2z == 'rip') then
 in0r = 0
 in0s = 0
 in0g = 0
 iliqskin = 0
!call dbzcalc(qv,rainmr,snowmr,graupelmr,tmp,p,refl, &
!             ly,lx,nz,in0r,in0s,in0g,iliqskin)

 do j=1,ly
 do i=1,lx
!  Compute the max value in the column
   max_refl(i,j)=maxval(refl(i,j,:))
 enddo ! i
 enddo ! j

else
 do j=1,ly
 do i=1,lx
   do k=1,nz

! Compute the basic reflectivity using RAMS reflectivity algorithm.

      if(c_m2z == 'rams') then

      	w=264083.11*(      rainmr(i,j,k)  &
          	     +0.2*(snowmr(i,j,k))  &
        	     +2.0* graupelmr(i,j,k)                )
       	w=max(1.,w)
      	refl(i,j,k)=17.8*alog10(w)

        if (refl(i,j,k).eq.0.) then
            refl(i,j,k)=-10.0
        endif 

! Compute the basic reflectivity using Kessler reflectivity algorithm.

      elseif (c_m2z == 'kessler') then

        refl(i,j,k) =17300.0 * &
                    (rho(i,j,k) * 1000.0 * &
                    MAX(0.0,rainmr(i,j,k)))**svnfrth

        ! Add the ice component
        refl(i,j,k)=refl(i,j,k) + &
                    38000.0*(rho(i,j,k) * 1000.0 * &
!                   MAX(0.0,icemr(i,j,k)+snowmr(i,j,k)+graupelmr(i,j,k)))**2.2
                    MAX(0.0,             snowmr(i,j,k)+graupelmr(i,j,k)))**2.2

        ! Convert to dBZ (with a lower limit of zero, given the max application)
        refl(i,j,k) = 10.*ALOG10(MAX(refl(i,j,k),1.0))

        ! Set zero reflectivity to -10
        if (refl(i,j,k).eq.0.) then
            refl(i,j,k)=-10.0
        endif 

! Computing the reflectivity by using the synpolrad code

      elseif (c_m2z == 'synp') then

        call versuch(rainmr(i,j,k),snowmr(i,j,k),graupelmr(i,j,k),tmp(i,j,k) &   ! I
                    ,zhhx(i,j,k)                                             &   ! O
                    ,ahhx(i,j,k)                                             &   ! O
                    ,zhvx(i,j,k),zvhx(i,j,k),zvvx(i,j,k)                     &   ! O
                    ,p(i,j,k)                                                &   ! I
                    ,qv(i,j,k)                                               &   ! I
                    ,delahvx(i,j,k),kkdp(i,j,k),rhohvx(i,j,k))                   ! O 

!       qr = specific rain   (kg/kg)
!       qs = specific snow (kg/kg)
!       qg = specific graupel (kg/kg)
!       t  = temperature (K)
!       p  = pressure (Pa) - druck
!       qv = specific water vapor (kg/kg) - qd

!       intent(out):
!       zhh = reflectivity-hh (Z), not in dBZ
!       zhv
!       zvh
!       zvv
!       kappa = extinction coeff (km-1) - ahhx
!       delahv   = specific differential attenuationa (dB km-1)
!       kkdp   = specific differential phase (degree km-1)
!       rhohv = hv-vh cross correlation coeffcient (n/a)

!       Convert to dBZ
        zvvxmax = MAX(zvvx(i,j,k),1.0)
        zhhxmax = MAX(zhhx(i,j,k),1.0)
        refl(i,j,k) = 10.*ALOG10(zhhxmax)
        zdr(i,j,k)  = 10.*ALOG10(zhhx(i,j,k)/zvvxmax)
        ldr(i,j,k)  = 10.*ALOG10(zvhx(i,j,k)/zhhxmax)

      endif ! c_m2z (reflectivity algorithm)

! Since we are going from the ground up, we can 
! check threshold and set echo top.

      if (refl(i,j,k) >= max_top_thresh) echo_tops(i,j)=hgt(i,j,k) 
   enddo

! Compute the max value in the column

   max_refl(i,j)=maxval(refl(i,j,:))
 enddo
 enddo

endif

! Convert microphysical precipitation mixing ratios to concentrations - multiplying by rho
rainmr(:,:,:)    = rainmr(:,:,:) * rho(:,:,:)
snowmr(:,:,:)    = snowmr(:,:,:) * rho(:,:,:)
! icemr(:,:,:)     = icemr(:,:,:) * rho(:,:,:)
graupelmr(:,:,:) = graupelmr(:,:,:) * rho(:,:,:)

return
end

!===============================================================================

subroutine lfm_ua_pcptype(nx,ny,nz,rainmr_prs,snowmr_prs,graupelmr_prs  &
                         ,pcptype_prs)

implicit none

real, parameter :: nonecode = 0.
real, parameter :: raincode = 1.
real, parameter :: snowcode = 2.
real, parameter :: zraincode =3.
real, parameter :: sleetcode = 4.
real, parameter :: hailcode = 5.
real, parameter :: drizzlecode = 6.
real, parameter :: zdrizzlecode = 7.
real, parameter :: rainsnowcode = 8.
real, parameter :: rainicecode = 9.

integer :: nx,ny,nz,i,j,k
real, dimension(nx,ny,nz) :: rainmr_prs,snowmr_prs,graupelmr_prs,pcptype_prs

pcptype_prs=nonecode

do k=1,nz
do j=1,ny
do i=1,nx
   if (rainmr_prs(i,j,k) > 0.) then
      if (snowmr_prs(i,j,k) > 0.) then
         if (graupelmr_prs(i,j,k) > snowmr_prs(i,j,k)) then
            pcptype_prs(i,j,k)=rainicecode
         else
            pcptype_prs(i,j,k)=rainsnowcode
         endif
      else
         if (graupelmr_prs(i,j,k) > 0.) then
            pcptype_prs(i,j,k)=rainicecode
         else
            pcptype_prs(i,j,k)=raincode
         endif
      endif
   else
      if (snowmr_prs(i,j,k) > 0.) then
         if (graupelmr_prs(i,j,k) > snowmr_prs(i,j,k)) then
            pcptype_prs(i,j,k)=sleetcode
         else
            pcptype_prs(i,j,k)=snowcode
         endif
      else
         if (graupelmr_prs(i,j,k) >  0) then
            pcptype_prs(i,j,k)=sleetcode
         else
            pcptype_prs(i,j,k)=nonecode
         endif
      endif
   endif
enddo
enddo
enddo

return
end

!===============================================================================

subroutine lfm_sfc_pcptype(nx,ny,nz,rainmr_sig,snowmr_sig,graupelmr_sig  &
                          ,pcptype_sfc)

implicit none

real, parameter :: nonecode = 0.
real, parameter :: raincode = 1.
real, parameter :: snowcode = 2.
real, parameter :: zraincode =3.
real, parameter :: sleetcode = 4.
real, parameter :: hailcode = 5.
real, parameter :: drizzlecode = 6.
real, parameter :: zdrizzlecode = 7.
real, parameter :: rainsnowcode = 8.
real, parameter :: rainicecode = 9.

integer :: nx,ny,nz,i,j,k
real, dimension(nx,ny) :: pcptype_sfc
real, dimension(nx,ny,nz) :: rainmr_sig,snowmr_sig,graupelmr_sig

pcptype_sfc=nonecode

do j=1,ny
do i=1,nx
   if (rainmr_sig(i,j,1) > 0.) then
      if (snowmr_sig(i,j,1) > 0.) then
         if (graupelmr_sig(i,j,1) > snowmr_sig(i,j,1)) then
            pcptype_sfc(i,j)=rainicecode
         else
            pcptype_sfc(i,j)=rainsnowcode
         endif
      else
         if (graupelmr_sig(i,j,1) > 0.) then
            pcptype_sfc(i,j)=rainicecode
         else
            pcptype_sfc(i,j)=raincode
         endif
      endif
   else
      if (snowmr_sig(i,j,1) > 0.) then
         if (graupelmr_sig(i,j,1) > snowmr_sig(i,j,1)) then
            pcptype_sfc(i,j)=sleetcode
         else
            pcptype_sfc(i,j)=snowcode
         endif
      else
         if (graupelmr_sig(i,j,1) > 0.) then
            pcptype_sfc(i,j)=sleetcode
         else
            pcptype_sfc(i,j)=nonecode
         endif
      endif
   endif
enddo
enddo

return
end

!===============================================================================

subroutine interp_press_to_z(press_levels,heights,new_level,press_out  &
                            ,nx,ny,nz)

! Given a 1D array of pressure levels and a 3D array of heights (m) at
! those levels, this routine interpolates the pressure to the desired
! new_level.
! Pressures are in Pa, heights are in m!

implicit none

integer :: nx,ny,nz,i,j,k

real :: press_levels(nz),heights(nx,ny,nz),new_level,press_out(nx,ny)  &
       ,ptop,pbot,ztop,zbot

do j=1,ny
do i=1,nx
do k=2,nz
   if (heights(i,j,k) >= new_level) then
      ztop=heights(i,j,k)
      zbot=heights(i,j,k-1)
      ptop=press_levels(k)*0.01    ! Convert to mb
      pbot=press_levels(k-1)*0.01
      press_out(i,j)=exp((new_level*alog(pbot/ptop)  &
                         -ztop*alog(pbot)            &
                         +zbot*alog(ptop))           &
                        /(zbot-ztop))*100.  ! convert back to Pa
      exit
   endif
enddo
enddo
enddo

return
end

!===============================================================================

subroutine helicity(usig,vsig,zsig,usf,vsf,terrain,imax,jmax,ksig,srelhel)

! NAME: HELICITY
!
! PURPOSE
! =======
! CALCULATE STORM RELATIVE HELICITY.
!
! REMARKS
! =======
! HELICITY IS EQUAL TO THE VERTICAL INTEGRAL OF:
!
! (V - c) x dV/dz
!
! WHERE V IS THE MODEL PREDICTED WIND AND C IS THE 
! ESTIMATED STORM MOTION VECTOR. 
!
! UPDATES
! =======
! ??? 97   INITIAL VERSION......................................DNXM
! DEC 98   MODIFIED CALL TO WDIR FUNCTION.  AN ERROR IN WDIR WAS 
!          CORRECTED............................................DNXM
! JAN 99   CHANGED VARIABLE NAMES USIGCRS AND VSIGCRS TO USIG AND
!          AND VSIG AS WINDS ARE NOW AT THE CROSS POINT.........DNXM
! JAN 01   Adapted for use at FSL by B. Shaw, NOAA/FSL

use constants

implicit none

integer :: i,imax,j,jmax,k,k3km,k10km,ksig
real :: dz,stormdir,stormspd,stormu,stormv,sumdz,sumhel,sumu,sumv,wdir,wspd,dudz,dvdz,udif,vdif,zdif
real, dimension(imax,jmax) :: srelhel,terrain,usf,vsf
real, dimension(imax,jmax,ksig) :: usig,vsig,zsig

! Determine indices of 3 and 10 km layers (above ground level).

do j=1,jmax
do i=1,imax
   do k=1,ksig

      if (((zsig(i,j,k)-terrain(i,j)) >= 3000.)) then
         k3km=k
         exit
      endif
   enddo

   do k=ksig,k3km,-1
      if (((zsig(i,j,k)-terrain(i,j)) <= 10000.)) then
         k10km=k
         exit 
      endif
   enddo

! Estimate storm motion vector.  Speed is 75 percent of the mean
!   wind between 3 and 10 km.  Direction is 30 degrees to the
!   right of the mean wind.

   sumu=0.
   sumv=0.
   sumdz=0.

   do k=k3km,k10km
      dz=zsig(i,j,k)-zsig(i,j,k-1)
      sumdz=sumdz+dz
      sumu=sumu+(0.5*dz*(usig(i,j,k)+usig(i,j,k-1)))
      sumv=sumv+(0.5*dz*(vsig(i,j,k)+vsig(i,j,k-1)))
   enddo
   stormu=sumu/sumdz
   stormv=sumv/sumdz
   stormspd=wspd(stormu,stormv)*0.75

! When calling wdir, send in a cone factor of zero
!   so that direction is grid relative.

   stormdir=wdir(stormu,stormv,0.,0.,0.)+30.
   if (stormdir > 360.) stormdir=stormdir-360.

   stormu=-stormspd*sin(stormdir*deg2rad)
   stormv=-stormspd*cos(stormdir*deg2rad)

   sumhel=0.
  
! Calculate helicity.  Integrate between the ground and 3 km,
!   A depth that is frequently used in the literature.

   do k=1,k3km
      if (k .eq. 1) then
       udif = 0.5*(usig(i,j,k)+usf(i,j))-stormu
       vdif = 0.5*(vsig(i,j,k)+vsf(i,j))-stormv
       dz=zsig(i,j,k)-terrain(i,j)
       dudz = (usig(i,j,k)-usf(i,j))/dz 
       dvdz = (vsig(i,j,k)-vsf(i,j))/dz 
      else
       udif = 0.5*(usig(i,j,k)+usig(i,j,k-1))-stormu
       vdif = 0.5*(vsig(i,j,k)+vsig(i,j,k-1))-stormv
       dz=zsig(i,j,k)-zsig(i,j,k-1)
       dudz = (usig(i,j,k)-usig(i,j,k-1))/dz 
       dvdz = (vsig(i,j,k)-vsig(i,j,k-1))/dz 
      endif
      sumhel=sumhel + (udif*dvdz - vdif*dudz)*dz
   enddo
   srelhel(i,j)=-sumhel
enddo
enddo

return
end

!===============================================================================

subroutine updraft_helicity(usig,vsig,wsig,zsig,zsigf,terrain,lat,lon,imax,jmax,ksig,uhel)

! updraft_helicity
!
! Calculate the helicity of updrafts within the 2-5km AGL layer of the updraft
!
! Updraft helicity is equal to the vertical integral bounded by 2-5km of
! W x zeta
! where W is the vertical velocity and zeta is the relative vorticity
! A centered-difference scheme is used.
!
! UPDATES
! =======
! 09 07   INITIAL VERSION.........................................CJA

use constants

implicit none

integer :: i,imax,j,jmax,k,k2km,k5km,ksig
real :: dz,dx,dy,dvdx,dudy,sumhel
real, dimension(imax,jmax) :: uhel,terrain,lat,lon
real, dimension(imax,jmax,ksig) :: wsig,usig,vsig,zsig
real, dimension(imax,jmax,ksig+1) :: zsigf

! Determine indices of 2 and 5 km layers (above ground level).

do j=2,jmax-1
do i=2,imax-1
   do k=1,ksig

      if (((zsig(i,j,k)-terrain(i,j)) >= 2000.)) then
         k2km=k
         exit
      endif
   enddo

   do k=ksig,k2km,-1
      if (((zsig(i,j,k)-terrain(i,j)) <= 5000.)) then
         k5km=k
         exit 
      endif
   enddo

! Calculate updraft helicity.  Integrate between 2 and 5 km.

   sumhel=0.
   do k=k2km,k5km
      dz = zsigf(i,j,k)-zsigf(i,j,k-1)
      dx = (lon(i+1,j)-lon(i-1,j))*cos(lat(i,j)*deg2rad)*eradius*deg2rad
      dvdx = (vsig(i+1,j,k)-vsig(i-1,j,k))/dx
      dy = (lat(i,j+1)-lat(i,j-1))*eradius*deg2rad
      dudy = (usig(i,j+1,k)-usig(i,j-1,k))/dy
      sumhel = sumhel + wsig(i,j,k)*(dvdx-dudy)*dz 
      if (k==k5km) then
!      print *, i,j,k,dx,dy
!      print *, dz,dvdx,dudy,wsig(i,j,k),sumhel
      endif
   enddo
   uhel(i,j)=sumhel
enddo
enddo

return
end


!===============================================================================

subroutine bulk_shear(usig,vsig,zsig,zsigf,terrain,lat,lon,imax,jmax,ksig,bshr)

! bulk_shear - sfc to 6km agl

use constants

implicit none

integer :: i,imax,j,jmax,k,k0km,k6km,ksig
real    :: usigl,usigh,vsigl,vsigh,ushr,vshr
real, dimension(imax,jmax) :: bshr,terrain,lat,lon
real, dimension(imax,jmax,ksig) :: usig,vsig,zsig
real, dimension(imax,jmax,ksig+1) :: zsigf

! Determine indices of 0 and 6 km layers (above ground level).

do j=1,jmax
do i=1,imax
   do k=1,ksig
      if (((zsig(i,j,k)-terrain(i,j)) >=    0.)) then
         k0km=k
         exit
      endif
   enddo

   do k=ksig,k0km,-1
      if (((zsig(i,j,k)-terrain(i,j)) <= 6000.)) then
         k6km=k
         exit 
      endif
   enddo

! Calculate bulk shear.  Difference between 0 and 6 km.

   usigl = usig(i,j,k0km)
   usigh = usig(i,j,k6km)
   vsigl = vsig(i,j,k0km)
   vsigh = vsig(i,j,k6km)

   ushr = usigh - usigl
   vshr = vsigh - vsigl

   bshr(i,j)=sqrt(ushr**2 + vshr**2)

enddo
enddo

return
end

!===============================================================================

subroutine ventilation(usig,vsig,zsig,pblhgt,topo,nx,ny,nz  &
                      ,upbl,vpbl,vent_ind)

implicit none

integer :: nx,ny,nz,i,j,k,nbl
real :: usum,vsum,umean,vmean,spmean
real, dimension(nx,ny) :: pblhgt,topo,upbl,vpbl,vent_ind
real, dimension(nx,ny,nz) :: usig,vsig,zsig

do j=1,ny
do i=1,nx
   if (pblhgt(i,j) > 0.) then

! Compute mean wind within the PBL.

      nbl=0
      usum=0.
      vsum=0.
      umean=0.
      vmean=0.
      do k=1,nz
         if (zsig(i,j,k)-topo(i,j) <= pblhgt(i,j)) then
            nbl=nbl+1
            usum=usum+usig(i,j,k)
            vsum=vsum+vsig(i,j,k)
         else
            exit
         endif
      enddo
      if (nbl > 0) then
         umean=usum/float(nbl)
         vmean=vsum/float(nbl)

! Compute mean wind speed for this layer.

         spmean=sqrt(umean**2+vmean**2)

! Multiply mean speed by PBL depth to get index.

         vent_ind(i,j)=pblhgt(i,j)*spmean
         upbl(i,j)=umean
         vpbl(i,j)=vmean
      else

! PBL height is lower than the lowest model level...use
! lowest model wind.

         spmean=sqrt(usig(i,j,1)**2+vsig(i,j,1)**2)
         vent_ind(i,j)=pblhgt(i,j)*spmean
         upbl(i,j)=usig(i,j,1)
         vpbl(i,j)=vsig(i,j,1)
      endif
   else
      if (pblhgt(i,j) < 0.)  &
         print*,'Warning:  PBL Height <0 in ventilation index:',pblhgt(i,j)
      vent_ind(i,j)=0.
      upbl(i,j)=0.
      vpbl(i,j)=0.
   endif
enddo
enddo

return
end

!===============================================================================

subroutine wind_agl(usig,vsig,zsig,aglhgt,topo,lx,ly,nz  &
                   ,uagl,vagl,r_missing_data)

implicit none

integer :: lx,ly,nz,i,j,k,nbl
logical :: l_debug
real :: aglhgt,refhgt,usum,vsum,umean,vmean,spmean,r_missing_data,zlow,zhigh
real, dimension(lx,ly) :: topo,uagl,vagl
real, dimension(lx,ly,nz) :: usig,vsig,zsig

write(6,*)' Subroutine wind_agl...'

do j=1,ly
do i=1,lx

! Compute mean wind within the PBL.
      if(i .eq. lx/2 .AND. j .eq. ly/2)then
         l_debug = .true.
      else
         l_debug = .false.
      endif       

      nbl=0
      usum=0.
      vsum=0.
      umean=0.
      vmean=0.
      refhgt = topo(i,j) + aglhgt
      do k=1,nz-1
         if (zsig(i,j,k) <= refhgt .AND. zsig(i,j,k+1) >= refhgt) then
            zlow = refhgt - zsig(i,j,k)
            zhigh = zsig(i,j,k+1) - refhgt
            nbl=nbl+1

            usum = 0.
            usum=usum+usig(i,j,k+1) * (zlow  / (zlow + zhigh))
            usum=usum+usig(i,j,k)   * (zhigh / (zlow + zhigh))

            vsum = 0.
            vsum=vsum+vsig(i,j,k+1) * (zlow  / (zlow + zhigh))
            vsum=vsum+vsig(i,j,k)   * (zhigh / (zlow + zhigh))

            exit
         endif
      enddo

      if (nbl > 0) then
         if(l_debug)then
             write(6,*)'topo,refhgt,zsig1,zsig2',topo(i,j),refhgt,zsig(i,j,k),zsig(i,j,k+1)
             write(6,*)'zlow,zhigh,usig(i,j,k),usig(i,j,kp1),usum = ',zlow,zhigh,usig(i,j,k),usig(i,j,k+1),usum
         endif
         uagl(i,j)=usum
         vagl(i,j)=vsum
      else
         uagl(i,j)=r_missing_data
         vagl(i,j)=r_missing_data
      endif

enddo
enddo

write(6,*)aglhgt,' meter wind U range = ',minval(uagl),maxval(uagl)
write(6,*)aglhgt,' meter wind V range = ',minval(vagl),maxval(vagl)

return
end

!===============================================================================

subroutine haines_layer(p3d,t3d_k,td3d_k,haines2d,nx,ny,nz  &
                       ,pmbbot,pmbtop)

! Computes haines index for layer bounded by top and bottom pressure 
!  levels pmbtop and pmbbot (mb).

implicit none

integer :: nx,ny,nz,i,j,k,kk,km1
real :: p3d(nx,ny,nz)     ! 3D pressure in Pa
real :: t3d_k (nx,ny,nz)  ! 3D Temp in K
real :: td3d_k(nx,ny,nz)  ! 3D Dewpoint in K
real :: haines2d(nx,ny)   ! 2D haines index
real :: pmbbot,pmbtop     ! Bounding mb levels
real :: ppabot,ppatop,factor
real :: tmkt,tmkb,tdkb,deltat,dpdep,factor1,factor2

ppabot=pmbbot*100.
ppatop=pmbtop*100.

do j=1,ny
do i=1,nx
   if (p3d(i,j,1) >= ppabot) then
      do k=1,nz-1
        
! Find temperature at the top.

         if (p3d(i,j,k) > ppatop .and. p3d(i,j,k+1) <= ppatop) then
            tmkt=t3d_k(i,j,k)+(t3d_k(i,j,k+1)-t3d_k(i,j,k))   &
                *(alog(ppatop/p3d(i,j,k)))                    &
                /(alog(p3d(i,j,k+1)/p3d(i,j,k)))
         endif

! Find Temp/dewpoint at the bottom of the layer.
            
         if (p3d(i,j,k) > ppabot .and. p3d(i,j,k+1) <= ppabot) then
            factor=(alog(ppabot      /p3d(i,j,k)))  &
                  /(alog(p3d(i,j,k+1)/p3d(i,j,k)))
            tmkb= t3d_k(i,j,k)+( t3d_k(i,j,k+1)- t3d_k(i,j,k))*factor
            tdkb=td3d_k(i,j,k)+(td3d_k(i,j,k+1)-td3d_k(i,j,k))*factor
         endif

      enddo

      deltat=tmkb-tmkt
      dpdep =tmkb-tdkb

      if (nint(pmbbot) == 700) then  ! High Haines
         if (deltat <= 17.5) then
            factor1=1.
         elseif (deltat > 17.5 .and. deltat <= 21.5 ) then
            factor1=2.   
         else
            factor1=3.
         endif

         if (dpdep <= 14.5) then
            factor2=1.
         elseif (dpdep > 14.5 .and. dpdep <= 20.5) then
            factor2=2.
         else  
            factor2=3.
         endif

      elseif (nint(pmbbot) == 850) then  ! Mid-level Haines
         if (deltat <= 5.5) then
            factor1=1.
         elseif (deltat > 5.5 .and. deltat <= 10.5) then
            factor1=2.
         else
            factor1=3.
         endif

         if (dpdep <= 5.5) then
            factor2=1.
         elseif (dpdep > 5.5 .and. dpdep <= 12.5) then
            factor2=2.
         else
            factor2=3.
         endif
   
      else 
         print*,'Haines_layer needs 850 or 700 as bottom layer'
         print*,'Bottom level (mb) specified:',pmbbot
         stop 'bad_haines_layer'

      endif

      haines2d(i,j)=factor1+factor2

   endif 

enddo
enddo

return
end

!===============================================================================

subroutine fosberg_fwi(t2k,rh2pct,u10,v10,nx,ny,fwi)

! Computes the Fosberg Fire Weather Index
 
implicit none

integer :: nx,ny,i,j
real :: uuu10,vvv10,m,n,t2f,rh2
real :: t2k(nx,ny)      ! Sfc Temp (K)
real :: rh2pct(nx,ny)   ! Sfc RH (%)
real :: u10(nx,ny)      ! 10 m U wind (m/s)
real :: v10(nx,ny)      ! 10 m V wind (m/s)
real :: fwi(nx,ny)      ! Fosberg Index

do j=1,ny
do i=1,nx

! Convert Temperature from K to F.

   t2f=1.8*(t2k(i,j)-273.15)+32.0
        
! Convert u/v from m/s to mph.

   uuu10=u10(i,j)*2.237
   vvv10=v10(i,j)*2.237
   rh2=rh2pct(i,j)

   if (rh2 <= 10.5) then
      m=0.03229+(0.281073*rh2)-(0.000578*rh2*t2f)
   elseif (rh2 > 10.5 .and. rh2 <= 50.5) then
      m=2.22749+(0.160107*rh2)-(0.014784*t2f)
   elseif (rh2 > 50.5 .and. rh2 <= 100.) then
      m=21.0606+(0.005565*rh2**2)-(0.00035*rh2*t2f)-(0.483199*rh2)
   else 
      m=21.0606+(0.005565*100**2)-(0.00035*100*t2f)-(0.483199*100)
   endif

   n=1.-2.*(m/30.)+1.5*(m/30.)**2-0.5*(m/30.)**3
   fwi(i,j)=(n*sqrt(1.+uuu10**2+vvv10**2))/0.3002

enddo
enddo

return
end

!===============================================================================

subroutine capecin(psig,tsig,thetaesig,thetasig                &
                  ,rhsig,zsigfull,tprs,li,posbuoyen,negbuoyen  &
                  ,k500,imax,jmax,ksig,kprs)

! Purpose
! =======
! Calculate positive buoyant energy (or convective available
! energy) and negative buoyant energy (or convective inhibition).
!
! Also calculate lifted index.
!
! Reference
! =========
! Doswell and Rasmussen (1994), Wea and Fcsting, p 625.
!
! Updates
! =======
! 4 Jan 01  - Adapted from USAF Weather Agency routine
!             B. Shaw, NOAA/FSL

use constants
use lfmgrid, ONLY: large_pgrid

implicit none

integer :: i,imax,j,jmax,k,k500,kmax,ksig,kprs
real, parameter :: pref=1000.
real :: prslcl,thw,thwmax,tlcl,tmplcl,tparcel,wobf,deltaz,dtheta,thetaparcel,potential_temp,psave
real, dimension(imax,jmax) :: li,negbuoyen,posbuoyen
real, dimension(imax,jmax,ksig) :: psig,rhsig,thetasig,thetaesig,tsig
real, dimension(imax,jmax,kprs) :: tprs
real, dimension(imax,jmax,ksig+1) :: zsigfull
real, allocatable, dimension(:) :: buoy
logical :: compute_cin

posbuoyen=0.
negbuoyen=0.
allocate(buoy(ksig)) 
do j=1,jmax
do i=1,imax
   thwmax=-9999.0

   do k=1,ksig

! Pick the most unstable parcel in the lowest 50 mb as
!   indicated by the sigma level with the highest wet bulb
!   potential temperature.  Store index in kmax.

      if (((psig(i,j,1)-psig(i,j,k)) < 50.)) then
         thw=thetaesig(i,j,k)-wobf(thetaesig(i,j,k)-t0)
         if (thw > thwmax) then
            kmax=k
            thwmax=thw
         endif
      else
         exit
      endif

   enddo

! Calculate lifted index by lifting the most unstable parcel.

   if(.not. large_pgrid)then
       call the2t(thetaesig(i,j,kmax),500.,tparcel)
       li(i,j)=tprs(i,j,k500)-tparcel
   endif
    
! Calculate the temperature and pressure of the lifting 
!   condensation level.

   tmplcl=tlcl(tsig(i,j,kmax),rhsig(i,j,kmax))
   prslcl=psig(i,j,kmax)*(tmplcl/tsig(i,j,kmax))**cpor

! Calculate the buoyancy.

   posbuoyen(i,j)=0.
   negbuoyen(i,j)=0.
   do k=kmax,ksig

! Above the lcl, calculate virtual temperature of the
!   parcel as it moves along a moist adiabat.  Below the
!   lcl, lift parcel along a dry adiabat.
            
      if (psig(i,j,k) <= prslcl) then
         call the2t(thetaesig(i,j,kmax),psig(i,j,k),tparcel)
      else
         tparcel=thetasig(i,j,kmax)/(pref/psig(i,j,k))**kappa
      endif

! Compute the potential temperature of the parcel.

      thetaparcel=potential_temp(tparcel,psig(i,j,k)*100.)
      dtheta=thetaparcel-thetasig(i,j,k)
      deltaz=zsigfull(i,j,k+1)-zsigfull(i,j,k)
      buoy(k)=deltaz*dtheta/thetasig(i,j,k)
   enddo
          
! Now loop through the column again, partitioning the buoyency
!   into positive (CAPE) and negative (CIN) component.  We terminate
!   the contribution to CIN when/if a layer of CAPE greater than
!   150 mb deep is found.
 
   compute_cin=.true.
   psave=-100.
   do k=kmax,ksig
      if (buoy(k) > 0.) then
         if (psave < 0) then
            psave=psig(i,j,k)
         else
            if ((psave-psig(i,j,k)) > 150.) then
               compute_cin = .false.
            endif
         endif
         posbuoyen(i,j)=posbuoyen(i,j)+buoy(k)
      elseif (buoy(k) < 0.) then
         psave=-100.
         if (compute_cin) THEN
             negbuoyen(i,j)=negbuoyen(i,j)+buoy(k)
         endif
      endif 
   enddo
enddo
enddo
    
posbuoyen=grav*posbuoyen
negbuoyen=grav*negbuoyen
 
! Cap the negative buoyancy to a maximum value of 700 J/kg

where(negbuoyen < -700) negbuoyen=-700. 
deallocate(buoy)

return
end

!===============================================================================

subroutine wintprec(tsig,zsig,zprs,psfc,tsfc,tdsfc,terrain,prcpinc,imax,jmax  &
                   ,ksig,kprs,k700,k850,k1000,prcpconv,preciptype)         

! (Winter) Precipitation Algorithm
!
! Purpose:  This product identifies areas of precipitation and the
! type expected, based on MM5 data.  The process is essentially two-
! fold.  First, areas of precipitation are identified when the MM5 
! precipitation array (PRCPINC) exceeds 0.01 inch. 
!
! Second, thickness thresholds are used at the gridpoints for which 
! it has been determined that precipitation is occurring and the    
! surface pressure is greater than 850mb (i.e., non-mountainous    
! regions).  The thickness thresholds utilized are based on         
! meteorological research from the pertinent sources
! (e.g., MWR, W&F), and are as follows: 
!
!                             (THICK1)     (THICK2)   
!                           1000mb-850mb  850mb-700mb  <--THICKNESS
!
! PRECIPITATION TYPE:   RAIN   GT 1310      GT 1540 
!
!              FREEZING RAIN   GT 1310      GT 1540 [sig1 T < 0Co]
!
!                  ICE/MIXED   LE 1310      GT 1540 
!
!                       SNOW   LE 1310      LE 1540. 
!
! Over mountainous terrain, precipitation type is limited to either 
! rain or snow.  This is consistent with climatic data presented in 
! "A Regional Climatology of Freezing Precipitation for the        
! Contiguous United States" (10th Conf. on Applied Climatology, 20-
! 24 Oct 97).  Where a precipitation occurrence has been determined, 
! the temperatures in the lowest 1500 m are checked.  If all are 
! below freezing, SNOW is forecasted; otherwise RAIN is predicted.
!
! Modification:  Added ability to predict regions where thunderstorm
! activity may occur.  Prior to exiting the main loop, a check is 
! made:  where rain is predicted AND the convective component of the
! precip exceeds 0.01", forecast for thunderstorms.  
!
! UPDATES
! =======
! JAN 2001 1998  INITIAL VERSION, ADAPTED FROM USAF WEATHER AGENCY
!     Brent Shaw, NOAA/Forecast Systems Lab

use constants

implicit none

integer :: i,imax,j,jmax,k,k1000,k850,k700,kprs,ksig,k1500
real :: tsfcf,thickhigh,thicklow,tsig1,tsig2,tsig3,fahren,td_c,tw_c,pres_mb,twet_fast
real, dimension(imax,jmax) :: prcpconv,prcpinc,preciptype,psfc,terrain,tsfc,tdsfc
real, dimension(imax,jmax,ksig) :: tsig,zsig
real, dimension(imax,jmax,kprs) :: zprs

do j=1,jmax
do i=1,imax

! The threshold for calculating precip type is 0.0001 meter
!  per time period.

   if (prcpinc(i,j) <= 0.0001) then
      preciptype(i,j)=0
   else

! Check the surface pressure to determine whether high or low-elevation 
!  logic is used.  850 mb is the cutoff.

! Low elevation grid point.

      if (psfc(i,j) > 85000.) then

! Calculate thicknesses that will be used to determine precip type.

         thicklow=zprs(i,j,k850)-zprs(i,j,k1000)
         thickhigh=zprs(i,j,k700)-zprs(i,j,k850)

! Rain, or if surface temperature is below freezing, freezing rain.  

         if ((thicklow > 1310.) .and. (thickhigh > 1540.)) then
            if (tsfc(i,j) >= t0) then  
               preciptype(i,j)=1  
            else
               preciptype(i,j)=3  
            endif

! Ice/mixed.

         elseif (thicklow <= 1310. .and. thickhigh > 1540.) then
            preciptype(i,j)=4  

! Rain or snow.

         elseif (thicklow <= 1310. .and. thickhigh <= 1540.) then
            tsfcf=fahren(tsfc(i,j)-t0) 
            td_c=tdsfc(i,j)-t0
            pres_mb=psfc(i,j) / 100.
            tw_c=twet_fast(tsfc(i,j)-t0,td_c,pres_mb)
!           if (tsfcf >= 37.) then
            if (tw_c >= +1.3) then
               preciptype(i,j)=1    
            else                                  
               preciptype(i,j)=5   
            endif

! Rain.

         else
            preciptype(i,j)=1
         endif 

! High terrain grid point.

      else

! Find top of 1500 m agl layer.
              
         do k=1,ksig
            if (zsig(i,j,k)-terrain(i,j) >= 1500.) then 
               k1500=k
               exit
            endif
         enddo    

! If the model top is ever lowered, the above code
!  could fail on the top of a high mountain.

         k1500=max(k1500,1)
 
! Find temperature at the bottom, top and middle of
!  1500 m agl layer.  for middle layer, recycle variable k1500.

         tsig1=tsig(i,j,1)
         tsig3=tsig(i,j,k1500)
         k1500=max(1,nint(float(k1500)/2.))
         tsig2=tsig(i,j,k1500)

! Snow.
         td_c=tdsfc(i,j)-t0
         pres_mb=psfc(i,j) / 100.
         tw_c=twet_fast(tsfc(i,j)-t0,td_c,pres_mb)

 !       if (tsig1 < t0 .and. tsig2 < t0 .and. tsig3 .lt. t0) then 
         if (tw_c <= +1.3) then
            preciptype(i,j)=5  

! Rain & check for thunderstorms.

         else
            preciptype(i,j)=1    
         endif 
      endif
       
   endif

   if (preciptype(i,j) == 1 .and. prcpconv(i,j) >= 0.001)  &
      preciptype(i,j)=2         ! Thunderstorm

enddo
enddo

return
end

!===============================================================================

subroutine snowfall(tsig,tsfc,prcpinc,preciptype,imax,jmax,ksig,snowinc,snowtot) 

! Name: Snow Accumulation Algorithm
!
! Purpose
! =======
! This algorithm calculates incremental snow accumulation and
! total snow accumulation.
!
! Method
! ======
! - If precip type is snow
!    - calculate column max temperature on dot
!    - calculate incremental precip on dot
!    - Obtain snow_to_rain ratio from LAPS analysis subroutine
! - If precip type is not snow
!    - snow accumulation is zero
!
! Variables
! =========
! Name             Type      I/O     Definition
! ----             ----      ---     ----------
! imax             integer  input   grid dimension, i direction
! jmax             integer  input   grid dimension, j direction    
! prcpinc          real     input   3hr precip accum
!                                   array(2-D)
! preciptype       integer  input   0 - is no precip
!                                     - is rain
!                                   2 - is snow
!                                     - is ice/mixed
! snowinc          real     output  3hr snow accum
!                                       array(2-D)
! snowtot          real     output  total snow accum
!                                      array(2-D)
! tsfc             real     input   surface temp array (2-D)
! tsfcf            real     local   surface temp in F
!
! Updates
! =======
! 20 Feb 98  Initial version..................Capt. John Lewis/DNXT
! 10 Nov 98  Changed preciptype flag for snow 
!            from 4 to 5; this corresponds 
!            with changes made to the precip 
!            type algorithm (wintprec.f)...Capt David Beberwyk/DNXT   
!  5 Jan 99  Removed interpolations to dot grid as mmpost now
!            operates on the cross grid........................DNXM
!  4 Jan 01  Adapted by FSL for use with LAPS.. B. Shaw, NOAA/FSL
!    Jan 12  Use column max temperature and analysis library routine
!            to obtain snow_to_rain ratio (S. Albers)

use constants

implicit none

integer :: i,imax,j,jmax,ksig
real :: fahren,tsfcf,tcolmax,snow_to_rain_ratio
real, dimension(imax,jmax) :: prcpinc,preciptype,snowinc,snowtot,tsfc,tdsfc
real, dimension(imax,jmax,ksig) :: tsig
      
do j=1,jmax
do i=1,imax    
            
! Check if precipitation type is snow.

   if (preciptype(i,j) == 2) then

      if(.true.)then ! Liquid equivalent of snow depends on column max temp
          tcolmax = maxval(tsig(i,j,:))
          snowinc(i,j) = snow_to_rain_ratio(tcolmax) * prcpinc(i,j)

      else ! Liquid equivalent of snow depends on surface temperature.
          tsfcf=fahren(tsfc(i,j)-t0)
          if (tsfcf >= 10.0) then
             snowinc(i,j)=10.*prcpinc(i,j)
          else
             snowinc(i,j)=15.*prcpinc(i,j)
          endif
      endif

! If precip type is not snow then snow accum is zero.

   else
      snowinc(i,j)=0.
   endif
            
enddo
enddo

! Update snow total.

snowtot=snowtot+snowinc

return
end

!===============================================================================

subroutine height_tw(pr,ht,tp,td,spr,sht,stp,std,smr,slp,threshold,ztw0,lx,ly,nz)

implicit none

integer :: lx,ly,nz,i,j,k,ku,kl
real, parameter :: lapse=0.0065
real :: tw,twl,twu,rat,pr0,zbot,ztop,t0,rh0,td0,relhum,dewpt,threshold
real, dimension(lx,ly,nz) :: pr,ht,tp,td
real, dimension(lx,ly) :: spr,sht,stp,std,smr,slp,ztw0

ku = -999
kl = -999

do j=1,ly
do i=1,lx
   if(maxval(tp(i,j,:)) .lt. 273.15+threshold)then ! subfreezing sigma column
       goto 2
   endif
   if(kl .ne. -999)then ! use the previously found bracketing levels for efficiency
     twu=tw(tp(i,j,ku),td(i,j,ku),pr(i,j,ku))
     do k=kl,kl
      twl=tw(tp(i,j,k  ),td(i,j,k  ),pr(i,j,k  ))
      if (twu <= 273.15+threshold .and. twl > 273.15+threshold) then
         rat=(twl-(273.15+threshold))/(twl-twu)
         pr0=pr(i,j,k)+rat*(pr(i,j,k+1)-pr(i,j,k))
         if(pr0 .lt. 0.)then
             write(6,*)' ERROR in height_tw: pr0 < 0. ',pr0
             go to 1
         endif
         zbot=ht(i,j,k)
         ztop=ht(i,j,k+1)
         rat=alog(pr0/pr(i,j,k))/alog(pr(i,j,k+1)/pr(i,j,k))
         ztw0(i,j)=zbot+rat*(ztop-zbot)
         goto 1
      endif
     enddo
   endif
   twu=tw(tp(i,j,nz),td(i,j,nz),pr(i,j,nz))
   do k=nz-1,1,-1
      twl=tw(tp(i,j,k  ),td(i,j,k  ),pr(i,j,k  ))
      if (twu <= 273.15+threshold .and. twl > 273.15+threshold) then
         rat=(twl-(273.15+threshold))/(twl-twu)
         pr0=pr(i,j,k)+rat*(pr(i,j,k+1)-pr(i,j,k))
         if(pr0 .lt. 0.)then
             write(6,*)' ERROR in height_tw: pr0 < 0. ',pr0
             go to 1
         endif
         zbot=ht(i,j,k)
         ztop=ht(i,j,k+1)
         rat=alog(pr0/pr(i,j,k))/alog(pr(i,j,k+1)/pr(i,j,k))
         ztw0(i,j)=zbot+rat*(ztop-zbot)
         kl = k
         ku = k+1
         goto 1
      endif
      twu=twl
   enddo
2  continue ! check surface wet bulb
   twl=tw(stp(i,j),std(i,j),spr(i,j))
   if (twu <= 273.15+threshold .and. twl > 273.15+threshold) then
      rat=(twl-(273.15+threshold))/(twl-twu)
      pr0=spr(i,j)+rat*(pr(i,j,1)-spr(i,j))
      zbot=sht(i,j)
      ztop=ht(i,j,1)
      rat=alog(pr0/spr(i,j))/alog(pr(i,j,1)/spr(i,j))
      ztw0(i,j)=zbot+rat*(ztop-zbot)
   else
! Height of wet-bulb = 0C is below ground.
! Calculate wet-bulb temp at height = 0 and interpolate.
      twu=twl
      t0=stp(i,j)+lapse*sht(i,j)
      rh0=min(relhum(t0,smr,slp(i,j)),1.)
      td0=dewpt(t0,rh0)
      twl=tw(t0,td0,slp(i,j))
      if (twu <= 273.15+threshold .and. twl > 273.15+threshold) then
         rat=(twl-(273.15+threshold))/(twl-twu)
         pr0=slp(i,j)+rat*(spr(i,j)-slp(i,j))
         ztop=sht(i,j)
         rat=alog(pr0/slp(i,j))/alog(spr(i,j)/slp(i,j))
         ztw0(i,j)=rat*ztop
      else
         ztw0(i,j)=0.
      endif
   endif
1  continue
enddo
enddo

return
end

subroutine get_ghi_ratio(i4time1,i4time2,lat,lon,ni,nj,ghi_ratio)

! Ratio of GHI at time 2 to GHI at time 1

include 'trigd.inc'

real lat(ni,nj)
real lon(ni,nj)

real solalt1(ni,nj)
real solalt2(ni,nj)

real ghi_ratio(ni,nj)

call get_solalt_2d(lat,lon,i4time1,ni,nj,solalt1)
call get_solalt_2d(lat,lon,i4time2,ni,nj,solalt2)

write(6,*)' Range of solalt1 is ',minval(solalt1),maxval(solalt1)
write(6,*)' Range of solalt2 is ',minval(solalt2),maxval(solalt2)

do i = 1,ni
do j = 1,nj

    floor = max(.01 * (1.+solalt1(i,j)/3.),0.)
    zenfrac1 = max(sind(solalt1(i,j)),floor)

    floor = max(.01 * (1.+solalt2(i,j)/3.),0.)
    zenfrac2 = max(sind(solalt2(i,j)),floor)

    if(zenfrac1 .gt. 0)then
        ratio = zenfrac2 / zenfrac1
    else
        ratio = 1.
    endif

!   if(i .eq. ni/2)then
!       write(6,*)'solalt1/2,zenfrac1,zenfrac2,ratio',i,j,solalt1(i,j),solalt2(i,j),zenfrac1,zenfrac2,ratio
!   endif

    ghi_ratio(i,j) = min(max(ratio,.1),10.)

enddo ! j
enddo ! i

return
end

       
