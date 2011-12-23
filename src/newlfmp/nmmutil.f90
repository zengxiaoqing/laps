module nmmutil

implicit none

save

integer :: ncid

end module

!===============================================================================

subroutine get_nmm_dims(fname,nx,ny,nz)

use nmmutil

implicit none

include 'netcdf.inc'

integer :: nx,ny,nz,nid,icode
character(len=*) :: fname
logical :: there

! Open nmm file, and leave open for future use.

inquire(file=trim(fname),exist=there)
if (there) then
   icode=nf_open(trim(fname),nf_nowrite,ncid)
   if (ncid <= 0) then
      print*,'Could not open nmm file: ',trim(fname)
      stop
   endif
else
   print*,'Could not find nmm file: ',trim(fname)
   stop
endif

! Read nmm grid dimensions.

icode=nf_inq_dimid(ncid,'west_east',nid)
icode=nf_inq_dimlen(ncid,nid,nx)
icode=nf_inq_dimid(ncid,'south_north',nid)
icode=nf_inq_dimlen(ncid,nid,ny)
icode=nf_inq_dimid(ncid,'bottom_top',nid)
icode=nf_inq_dimlen(ncid,nid,nz)

return
end

!===============================================================================

subroutine fill_nmm_grid

use lfmgrid
use nmmutil
use constants

implicit none

include 'netcdf.inc'

integer :: nid,icode,k
real, allocatable, dimension(:,:) :: fld2d
real, allocatable, dimension(:,:,:) :: fld3d

nprojection='ROTATED LAT-LON'
icode=nf_get_att_real(ncid,NF_GLOBAL,'CEN_LAT',ntruelat1)
icode=nf_get_att_real(ncid,NF_GLOBAL,'CEN_LON',nstdlon)
icode=nf_get_att_real(ncid,NF_GLOBAL,'DX',ngrid_spacingx)
icode=nf_get_att_real(ncid,NF_GLOBAL,'DY',ngrid_spacingy)

icode=nf_inq_varid(ncid,'FIS',nid)
icode=nf_get_var_real(ncid,nid,nzsfc)
nzsfc(:,:)=nzsfc(:,:)/grav
icode=nf_inq_varid(ncid,'GLAT',nid)
icode=nf_get_var_real(ncid,nid,nlat)
icode=nf_inq_varid(ncid,'GLON',nid)
icode=nf_get_var_real(ncid,nid,nlon)

allocate(fld3d(nx,ny,nz+1))
icode=nf_inq_varid(ncid,'PINT',nid)
icode=nf_get_var_real(ncid,nid,fld3d)
icode=nf_inq_varid(ncid,'T',nid)
icode=nf_get_var_real(ncid,nid,ntsig)
icode=nf_inq_varid(ncid,'Q',nid)
icode=nf_get_var_real(ncid,nid,nmrsig)

! Compute height and pressure at half-levels.
call nmm_height(nx,ny,nz,fld3d,ntsig,nmrsig,nzsfc,npsig,nzsig)

icode=nf_inq_varid(ncid,'U',nid)
if (icode .eq. 0) icode=nf_get_var_real(ncid,nid,nusig)
icode=nf_inq_varid(ncid,'V',nid)
if (icode .eq. 0) icode=nf_get_var_real(ncid,nid,nvsig)
icode=nf_inq_varid(ncid,'W',nid)
if (icode .eq. 0) then
   icode=nf_get_var_real(ncid,nid,fld3d)
   do k=1,nz
      nwsig(:,:,k)=(fld3d(:,:,k)+fld3d(:,:,k+1))*0.5
   enddo
   nwsfc(:,:)=fld3d(:,:,1)
endif

deallocate(fld3d)

! NMM cloud water fields:
!   CWM is total condensate.
!   F_ICE is fraction of CWM that is ice.
!   liquid water content is CWM - F_ICE*CWM.
!   F_RAIN is fraction of liquid water content that is rain.
!   F_RIMEF is ratio of total ice to unrimed ice (>=1).
!   Snow cannot be differentiated from ice in nmm output.

icode=nf_inq_varid(ncid,'CWM',nid)
if (icode .eq. 0) then
   icode=nf_get_var_real(ncid,nid,ncldicemr_sig) ! CWM
   icode=nf_inq_varid(ncid,'F_ICE',nid)
   if (icode .eq. 0) then
      icode=nf_get_var_real(ncid,nid,nsnowmr_sig) ! F_ICE
      icode=nf_inq_varid(ncid,'F_RAIN',nid)
      if (icode .eq. 0) then
         icode=nf_get_var_real(ncid,nid,nrainmr_sig) ! F_RAIN
         icode=nf_inq_varid(ncid,'F_RIMEF',nid)
         if (icode .eq. 0) then
            icode=nf_get_var_real(ncid,nid,ngraupelmr_sig) ! F_RIMEF
!           ncldliqmr_sig=ncldicemr_sig-ncldicemr_sig*nsnowmr_sig
            ncldliqmr_sig=ncldicemr_sig * (1.0 - nsnowmr_sig)

            nrainmr_sig=ncldliqmr_sig*nrainmr_sig

            ncldicemr_sig=ncldicemr_sig*nsnowmr_sig

            ngraupelmr_sig=ncldicemr_sig/ngraupelmr_sig

            nsnowmr_sig=rmsg
         else
            ncldliqmr_sig=rmsg
            nrainmr_sig=rmsg
            ncldicemr_sig=rmsg
            ngraupelmr_sig=rmsg
            nsnowmr_sig=rmsg
         endif
      endif
   endif
endif

! Accumulate precip.

!icode=nf_inq_varid(ncid,'CUPREC',nid)
!if (icode .eq. 0) then
!   icode=nf_get_var_real(ncid,nid,npcp_tot)
!else
!   npcp_tot=0.
!endif

! In the nmm, ACPREC is the total precipitation.
! That is, ACPREC is the sum of cumulus and explicit precipitation.
npcp_tot=0.
allocate(fld2d(nx,ny))
icode=nf_inq_varid(ncid,'ACPREC',nid)
if (icode .eq. 0) then
   icode=nf_get_var_real(ncid,nid,fld2d)
   npcp_tot=npcp_tot+fld2d
endif
deallocate(fld2d)

icode=nf_inq_varid(ncid,'PSHLTR',nid)
if (icode .eq. 0) icode=nf_get_var_real(ncid,nid,npsfc)
icode=nf_inq_varid(ncid,'TSHLTR',nid)
if ((icode .eq. 0) .and. (minval(npsfc) < rmsg)) then
   icode=nf_get_var_real(ncid,nid,ntsfc)
   ntsfc=ntsfc*(npsfc/p0)**kappa
endif
icode=nf_inq_varid(ncid,'QSHLTR',nid)
if (icode .eq. 0) icode=nf_get_var_real(ncid,nid,nmrsfc)
icode=nf_inq_varid(ncid,'U10',nid)
if (icode .eq. 0) icode=nf_get_var_real(ncid,nid,nusfc)
icode=nf_inq_varid(ncid,'V10',nid)
if (icode .eq. 0) icode=nf_get_var_real(ncid,nid,nvsfc)
icode=nf_inq_varid(ncid,'TGROUND',nid)
if (icode .eq. 0) icode=nf_get_var_real(ncid,nid,nground_t)
icode=nf_inq_varid(ncid,'SFCSHX',nid)
if (icode .eq. 0) icode=nf_get_var_real(ncid,nid,nshflux)
icode=nf_inq_varid(ncid,'SFCLHX',nid)
if (icode .eq. 0) icode=nf_get_var_real(ncid,nid,nlhflux)
icode=nf_inq_varid(ncid,'PBLH',nid)
if (icode .eq. 0) icode=nf_get_var_real(ncid,nid,npblhgt)
icode=nf_inq_varid(ncid,'RLWTOA',nid)
if (icode .eq. 0) icode=nf_get_var_real(ncid,nid,nlwout)

if (fcsttime == 0) then
   nshflux=rmsg
   nlhflux=rmsg
   npblhgt=rmsg
endif
if (maxval(npblhgt) < 1.) npblhgt=rmsg

return
end

!===============================================================================

subroutine nmm_height(nx,ny,nz,npri,ntp,nmr,sht,npr,nht)

use constants

implicit none

!real, parameter :: r=287.,g=9.8,rog=r/g

integer :: nx,ny,nz,i,j,k

real :: alnp(nx,ny)
real, dimension(nx,ny) :: sht
real, dimension(nx,ny,nz) :: ntp,nmr,npr,nht
real, dimension(nx,ny,nz+1) :: npri,nhti

! Compute pressure at half levels in Pa.

do k=1,nz
   npr(:,:,k)=(npri(:,:,k)+npri(:,:,k+1))*0.5
enddo

! Calculate height at full levels using hypsometric eqn.

nhti(:,:,1)=sht(:,:)
do k=1,nz
   nhti(:,:,k+1)=nhti(:,:,k)+rog*(ntp(:,:,k)+0.61*nmr(:,:,k))  &
                            *alog(npri(:,:,k)/npri(:,:,k+1))
enddo

! Logarithmic interpolation to half-levels.

do k=1,nz
   alnp(:,:)=alog((npr(:,:,k)   /npri(:,:,k)))  &
            /alog (npri(:,:,k+1)/npri(:,:,k))
   nht(:,:,k)=nhti(:,:,k)+(nhti(:,:,k+1)-nhti(:,:,k))*alnp(:,:)
enddo

return
end
