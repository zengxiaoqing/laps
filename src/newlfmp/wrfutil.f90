module wrfutil

implicit none

save

integer :: ncid

end module

!===============================================================================

subroutine get_wrf_dims(fname,nx,ny,nz,istatus)

use wrfutil

implicit none

include 'netcdf.inc'

integer :: nx,ny,nz,nid,icode,nRec,istatus
character(len=*) :: fname
logical :: there

! istatus: 1=good return, 0=ERROR/bad return
istatus = 1  ! assume good return

! Open wrf file, and leave open for future use.

inquire(file=trim(fname),exist=there)
if (there) then
   icode=nf_open(trim(fname),nf_nowrite,ncid)
   if (ncid <= 0) then
      print*,'Could not open wrf file: ',trim(fname)
      istatus = 0
      return
   else
      print*,'Opened wrf file: ',trim(fname)
   endif
else
   print*,'Could not find wrf file: ',trim(fname)
   istatus = 0
   return
endif

icode=nf_inq_dimid(ncid,'Time',nid)
icode=nf_inq_dimlen(ncid,nid,nRec)
print*,'Time (number of records): ', nRec

! Verify that file has 1 or more records (dimension "Time")
if (nRec .eq. 0) then
   print*,'ERROR: wrf file contains no data: ',trim(fname)
   print*,'ERROR: cannot create output files...STOPPING!'
   istatus = 0
   return
endif

! Read wrf grid dimensions.

icode=nf_inq_dimid(ncid,'west_east',nid)
icode=nf_inq_dimlen(ncid,nid,nx)
      print*,'west_east: ',nx
icode=nf_inq_dimid(ncid,'south_north',nid)
icode=nf_inq_dimlen(ncid,nid,ny)
      print*,'south_north: ',ny
icode=nf_inq_dimid(ncid,'bottom_top',nid)
icode=nf_inq_dimlen(ncid,nid,nz)
      print*,'bottom_top: ',nz

return
end

!===============================================================================

subroutine fill_wrf_grid

use lfmgrid
use wrfutil
use constants

implicit none

include 'netcdf.inc'

integer :: nid,icode,mapproj,i,j,k
real, allocatable, dimension(:,:) :: ncon_pcp_tot
real, allocatable, dimension(:,:,:) :: fld3d

! Fill native map projection settings.

icode=nf_get_att_int(ncid,NF_GLOBAL,'MAP_PROJ',mapproj)
icode=nf_get_att_real(ncid,NF_GLOBAL,'DX',ngrid_spacingx)
icode=nf_get_att_real(ncid,NF_GLOBAL,'DY',ngrid_spacingy)
icode=nf_get_att_real(ncid,NF_GLOBAL,'TRUELAT1',ntruelat1)
icode=nf_get_att_real(ncid,NF_GLOBAL,'TRUELAT2',ntruelat2)
icode=nf_get_att_real(ncid,NF_GLOBAL,'STAND_LON',nstdlon)

select case (mapproj)
   case(1)
      nprojection='LAMBERT CONFORMAL'
   case(2)
      nprojection='POLAR STEREOGRAPHIC'
   case(3)
      nprojection='MERCATOR'
end select

! Allocate local variables.

allocate(ncon_pcp_tot(nx,ny))

! Read model data.

if(l_process_uv)then ! U,V
  print*,' Reading WRF U/V ',large_ngrid,large_pgrid,l_process_uv
  icode=nf_inq_varid(ncid,'U',nid)
  if (icode .eq. 0) then
   allocate(fld3d(nx+1,ny,nz))
   icode=nf_get_var_real(ncid,nid,fld3d)
      print*,'U: ',icode
   do i=1,nx-1
      nusig(i,:,:)=(fld3d(i,:,:)+fld3d(i+1,:,:))*0.5
   enddo
   do i=nx,nx
      nusig(i,:,:)=(fld3d(i,:,:)+fld3d(i  ,:,:))*0.5
   enddo
   deallocate(fld3d)
  endif

  icode=nf_inq_varid(ncid,'V',nid)
  if (icode .eq. 0) then
   allocate(fld3d(nx,ny+1,nz))
   icode=nf_get_var_real(ncid,nid,fld3d)
      print*,'V: ',icode
   do j=1,ny-1
      nvsig(:,j,:)=(fld3d(:,j,:)+fld3d(:,j+1,:))*0.5
   enddo
   do j=ny,ny
      nvsig(:,j,:)=(fld3d(:,j,:)+fld3d(:,j  ,:))*0.5
   enddo
   deallocate(fld3d)
  endif
else
  print*,' Skipping read of WRF U/V ',large_ngrid,large_pgrid,l_process_uv
endif

icode=nf_inq_varid(ncid,'P',nid)
if (icode .eq. 0) then
   icode=nf_get_var_real(ncid,nid,npsig)
   icode=nf_inq_varid(ncid,'PB',nid)
   if (icode .eq. 0) then
    allocate(fld3d(nx,ny,nz))
      icode=nf_get_var_real(ncid,nid,fld3d)
      npsig=npsig+fld3d
      deallocate(fld3d)
   else
      npsig=rmsg
   endif
endif

icode=nf_inq_varid(ncid,'T',nid)
if (icode .eq. 0) then
   icode=nf_get_var_real(ncid,nid,ntsig)
   ntsig=(ntsig+300.)*(npsig/p0)**kappa
endif

icode=nf_inq_varid(ncid,'QVAPOR',nid)
if (icode .eq. 0) icode=nf_get_var_real(ncid,nid,nmrsig)

icode=nf_inq_varid(ncid,'QCLOUD',nid)
if (icode .eq. 0) icode=nf_get_var_real(ncid,nid,ncldliqmr_sig)

icode=nf_inq_varid(ncid,'QICE',nid)
if (icode .eq. 0) icode=nf_get_var_real(ncid,nid,ncldicemr_sig)

icode=nf_inq_varid(ncid,'QRAIN',nid)
if (icode .eq. 0) icode=nf_get_var_real(ncid,nid,nrainmr_sig)

icode=nf_inq_varid(ncid,'QSNOW',nid)
if (icode .eq. 0) icode=nf_get_var_real(ncid,nid,nsnowmr_sig)

icode=nf_inq_varid(ncid,'QGRAUP',nid)
if (icode .eq. 0) icode=nf_get_var_real(ncid,nid,ngraupelmr_sig)

if(c_m2z .eq. 'wrf')then
  print*,' Reading WRF REFL_10CM'
  icode=nf_inq_varid(ncid,'REFL_10CM',nid)
  if (icode .eq. 0) icode=nf_get_var_real(ncid,nid,nrefl_sig)
endif

if(l_process_w)then ! W   
  print*,' Reading WRF W ',large_ngrid,large_pgrid,l_process_w
  allocate(fld3d(nx,ny,nz+1))
  icode=nf_inq_varid(ncid,'W',nid)
  if (icode .eq. 0) then
   icode=nf_get_var_real(ncid,nid,fld3d)
   do k=1,nz
      nwsig(:,:,k)=(fld3d(:,:,k)+fld3d(:,:,k+1))*0.5
   enddo
  endif
  deallocate(fld3d)
else
  print*,' Skipping read of WRF W ',large_ngrid,large_pgrid,l_process_w
endif ! l_process_w

if(.true.)then ! Z
  print*,' Reading WRF Z'
  allocate(fld3d(nx,ny,nz+1))
  icode=nf_inq_varid(ncid,'PH',nid)
  if (icode .eq. 0) then
   icode=nf_get_var_real(ncid,nid,fld3d)
   do k=1,nz
      nzsig(:,:,k)=(fld3d(:,:,k)+fld3d(:,:,k+1))*0.5
   enddo
   icode=nf_inq_varid(ncid,'PHB',nid)
   if (icode .eq. 0) then
      icode=nf_get_var_real(ncid,nid,fld3d)
      do k=1,nz
         nzsig(:,:,k)=(nzsig(:,:,k)+(fld3d(:,:,k)+fld3d(:,:,k+1))*0.5)/grav
      enddo
   else
      nzsig=rmsg
   endif
  endif
  deallocate(fld3d) 
else
  print*,' Skipping read of WRF Z'
endif ! large_ngrid

icode=nf_inq_varid(ncid,'TSK',nid)
if (icode .eq. 0) icode=nf_get_var_real(ncid,nid,nground_t)

icode=nf_inq_varid(ncid,'PSFC',nid)
if (icode .eq. 0) icode=nf_get_var_real(ncid,nid,npsfc)

icode=nf_inq_varid(ncid,'T2',nid)
if (icode .eq. 0) icode=nf_get_var_real(ncid,nid,ntsfc)

icode=nf_inq_varid(ncid,'Q2',nid)
if (icode .eq. 0) icode=nf_get_var_real(ncid,nid,nmrsfc)

icode=nf_inq_varid(ncid,'U10',nid)
if (icode .eq. 0) icode=nf_get_var_real(ncid,nid,nusfc)

icode=nf_inq_varid(ncid,'V10',nid)
if (icode .eq. 0) icode=nf_get_var_real(ncid,nid,nvsfc)

icode=nf_inq_varid(ncid,'RAINC',nid)
if (icode .eq. 0) then
   icode=nf_get_var_real(ncid,nid,ncon_pcp_tot)
else
   ncon_pcp_tot=0.
endif

icode=nf_inq_varid(ncid,'RAINNC',nid)
if (icode .eq. 0) then
   icode=nf_get_var_real(ncid,nid,npcp_tot)
else
   npcp_tot=0.
endif

icode=nf_inq_varid(ncid,'HGT',nid)
if (icode .eq. 0) icode=nf_get_var_real(ncid,nid,nzsfc)

icode=nf_inq_varid(ncid,'XLAT',nid)
if (icode .eq. 0) icode=nf_get_var_real(ncid,nid,nlat)
      print*,'XLAT: ',icode,maxval(nlat)

icode=nf_inq_varid(ncid,'XLONG',nid)
if (icode .eq. 0) icode=nf_get_var_real(ncid,nid,nlon)
      print*,'XLONG: ',icode,maxval(nlon)

icode=nf_inq_varid(ncid,'SWDOWN',nid)
if (icode .eq. 0) icode=nf_get_var_real(ncid,nid,nswdown)

icode=nf_inq_varid(ncid,'GLW',nid)
if (icode .eq. 0) icode=nf_get_var_real(ncid,nid,nlwdown)

icode=nf_inq_varid(ncid,'OLR',nid)
if (icode .eq. 0) icode=nf_get_var_real(ncid,nid,nlwout)

icode=nf_inq_varid(ncid,'LH',nid)
if (icode .eq. 0) icode=nf_get_var_real(ncid,nid,nlhflux)

icode=nf_inq_varid(ncid,'GRDFLX',nid)
if (icode .eq. 0) icode=nf_get_var_real(ncid,nid,nshflux)

icode=nf_inq_varid(ncid,'PBLH',nid)
if (icode .eq. 0) then
   icode=nf_get_var_real(ncid,nid,npblhgt)
   if (maxval(npblhgt) <= 1.) npblhgt=rmsg
endif

icode=nf_close(ncid)

! Fill total precip and convert from mm to m.

npcp_tot=(npcp_tot+ncon_pcp_tot)*0.001

deallocate(ncon_pcp_tot)

if (fcsttime == 0.) then
   nlwdown=rmsg
   nswdown=rmsg
   nshflux=rmsg
   nlhflux=rmsg
endif

return
end
