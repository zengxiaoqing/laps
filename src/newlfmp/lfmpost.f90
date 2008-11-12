program lfmpost

use lfmgrid
use grib

implicit none

integer, external :: iargc
integer :: narg,chr,nc,i,j,k
integer :: laps_reftime,laps_valtime,istatus

character(len=256) :: laps_data_root,lfmprd_dir
character(len=24) :: a24time
character(len=16) :: atime,afcst
character(len=9) :: a9time
character(len=4) :: fcst_hhmm
character(len=2) :: adomnum,domnum_str

!beka 
integer ISTAT, I4_elapsed, ishow_timer, init_timer

logical :: back

! Read command line arguments.

!beka

	ISTAT = INIT_TIMER()


narg=iargc()
if (narg /= 6) call usage

!call ugetarg(1,mtype)
!call ugetarg(2,filename)
!call ugetarg(3,adomnum)
!call ugetarg(4,atime)
!call ugetarg(5,afcst)
!call ugetarg(6,laps_data_root)
call getarg(1,mtype)
call getarg(2,filename)
call getarg(3,adomnum)
call getarg(4,atime)
call getarg(5,afcst)
call getarg(6,laps_data_root)
call right_justify(adomnum)
call right_justify(atime)
call right_justify(afcst)
read(adomnum,'(i2)',err=900) domnum
read(atime,'(i16)',err=900) laps_reftime
read(afcst,'(i16)',err=900) laps_valtime
fcsttime=laps_valtime
laps_valtime=laps_reftime+laps_valtime

! Convert mtype to lower case and check for valid model type.

do i=1,len_trim(mtype)
   chr=ichar(mtype(i:i))
   if (chr > 64 .and. chr < 91) mtype(i:i)=char(chr+32)
enddo

if (trim(mtype) /= 'mm5' .and. trim(mtype) /= 'wrf' .and.  &
    trim(mtype) /= 'nmm' .and. trim(mtype) /= 'st4') then
   print *, 'Unrecognized model type provided as first arg: ',trim(mtype)
   print *, '   mm5, wrf, nmm, and st4 are supported.'
   stop
endif

print*,' '
print*,'LAPS forecast model postprocessing for ',trim(mtype),' file:'
print*,'   '//trim(filename)

! Read namelist.

back=.true.
nc=index(filename,'/',back)
if (nc > 1) then
   lfmprd_dir=filename(1:nc-1)
else
   lfmprd_dir="."
endif

call lfm_namelist(lfmprd_dir)

! Obtain native model grid dimensions.

call get_native_dims(mtype,filename,nx,ny,nz)

! Obtain output laps isobaric grid dimensions, pressure levels, and lat/lons.

call get_laps_static(laps_data_root)

! Allocate and fill input native model data.

call alloc_native_grid
call fill_native_grid
if (verbose) then
   print*,' '
   print*,'Native map projection parameters          Output map projection parameters'
   print*,'--------------------------------          --------------------------------'
   print'(a,3i5,11x,a,3i5)',' grid dimensions:',nx,ny,nz,'grid dimensions:',lx,ly,lz
   print'(a,f10.1,16x,a,f10.1)',' grid spacing   :',ngrid_spacingx,'grid spacing   :',grid_spacing
   print'(1x,2(''proj: '',a32,4x))',nprojection,projection
   print'(a,f10.3,16x,a,f10.3)',' true latitude 1:',ntruelat1,'true latitude 1:',truelat1
   print'(a,f10.3,16x,a,f10.3)',' true latitude 2:',ntruelat2,'true latitude 2:',truelat2
   print'(a,f10.3,16x,a,f10.3)',' std longitude  :',nstdlon,'std longitude  :',stdlon
   print*,' '
   if (trim(mtype) /= 'st4') then
     print'(a,2f10.2)', ' Min/Max value of native terrain: ',minval(nzsfc),maxval(nzsfc)
   endif
   print*,' '
   print*,'Corner points from native grid:'
   print*,'==============================='
   print*,' '
   print'(f8.3,1x,f8.3,10x,f8.3,1x,f8.3)'  &
        ,nlat(1,ny),nlon(1,ny),nlat(nx,ny),nlon(nx,ny)
   print*,'      (NW)----------------------(NE)'
   print*,'        |                        |'
   print*,'        |                        |'
   print*,'      (SW)----------------------(SE)'
   print'(f8.3,1x,f8.3,10x,f8.3,1x,f8.3)'  &
        ,nlat(1,1),nlon(1,1),nlat(nx,1),nlon(nx,1)
   print*, ' '
   if (trim(mtype) /= 'st4') then
     print*,'Diagnostics from native domain center:'
     print'(a,f7.0)',' Terrain height at center = ',nzsfc(nx/2,ny/2)
     print*,'------------------------------------------------------------'
     print*,'LEVEL     PRES(Pa)  HEIGHT     T      QV         U       V'
     print*,'------------------------------------------------------------'
     do k=1,nz
       print '(i5,2x,f12.1,2x,f6.0,2x,f6.2,2x,f8.6,2f8.2)'   &
             ,k,npsig(nx/2,ny/2,k),nzsig(nx/2,ny/2,k)  &
             ,ntsig(nx/2,ny/2,k),nmrsig(nx/2,ny/2,k)   &
             ,nusig(nx/2,ny/2,k),nvsig(nx/2,ny/2,k)
    enddo
  endif
endif

! Allocate and horizontally interpolate data.

call alloc_surface_grid
if (trim(mtype) /= 'st4') call alloc_hinterp_grid
call lfm_hinterp

call dealloc_grid('native')

!cj Added 6/21/2007
!cj Convert mixing ratio to specific humidity
hmrsig = hmrsig/(1.+hmrsig)

if (verbose) then
   print*,' '
   print'(a,2f10.2)', ' Min/Max value of hinterp terrain: ', minval(zsfc),maxval(zsfc)
   print*,' '
   print*,'Corner points from hinterp grid:'
   print*,'==============================='
   print*,' '
   print'(f8.3,1x,f8.3,10x,f8.3,1x,f8.3)'  &
        ,llat(1,ly),llon(1,ly),llat(lx,ly),llon(lx,ly)
   print*,'      (NW)----------------------(NE)'
   print*,'        |                        |'
   print*,'        |                        |'
   print*,'      (SW)----------------------(SE)'
   print'(f8.3,1x,f8.3,10x,f8.3,1x,f8.3)'  &
        ,llat(1,1),llon(1,1),llat(lx,1),llon(lx,1)
   if (trim(mtype) /= 'st4') then
     print*, ' '
     print*,'Diagnostics from hinterp domain center:'
     print'(a,f7.0)',' Terrain height at center = ',zsfc(lx/2,ly/2)
     print*,'------------------------------------------------------------'
     print*,'LEVEL     PRES(Pa)  HEIGHT     T      QV         U       V'
     print*,'------------------------------------------------------------'
     do k=1,nz
        print '(i5,2x,f12.1,2x,f6.0,2x,f6.2,2x,f8.6,2f8.2)'   &
              ,k,hpsig(lx/2,ly/2,k),hzsig(lx/2,ly/2,k)  &
              ,htsig(nx/2,ly/2,k),hmrsig(nx/2,ly/2,k)   &
              ,husig(nx/2,ly/2,k),hvsig(nx/2,ly/2,k)
     enddo
   endif
endif

! Set to zero any missing microphysics fields.

!if (make_micro) then
!   if (minval(hcldliqmr_sig)  == rmsg) hcldliqmr_sig=0.
!   if (minval(hcldicemr_sig)  == rmsg) hcldicemr_sig=0.
!   if (minval(hrainmr_sig)    == rmsg) hrainmr_sig=0.
!   if (minval(hsnowmr_sig)    == rmsg) hsnowmr_sig=0.
!   if (minval(hgraupelmr_sig) == rmsg) hgraupelmr_sig=0.
!endif

! Allocate and vertically interpolate data to isobaric grid.

if (trim(mtype) /= 'st4') then
 call alloc_isobaric_grid

!beka
	I4_elapsed = ishow_timer()

 call lfm_vinterp

!beka
	I4_elapsed = ishow_timer()

 if (verbose) then
    print*,' '
    print*,'Diagnostics from isobaric domain center:'
    print*,'--------------------------------------------'
    print*,'LEVEL     PRES(Pa)  HEIGHT     T      QV'
    print*,'--------------------------------------------'
    do k=1,lz
       print '(i5,2x,f12.1,2x,f6.0,2x,f6.2,2x,f8.6)'   &
             ,k,lprs(k)*100.,zprs(lx/2,ly/2,k)  &
             ,tprs(lx/2,ly/2,k),shprs(lx/2,ly/2,k)
    enddo
 endif

! Fill any missing surface fields and generate derived fields.

 call lfm_derived
 call lfm_vinterp
 if (verbose) then
    print*,' '
    print*,'lfm_reflectivity pressure max/min: ',maxval(refl_prs),minval(refl_prs)
    print*,'lfm_reflectivity hsig max/min: ',maxval(hrefl_sig),minval(hrefl_sig)
 endif

! Create point forecasts, if requested.

 call set_laps_projection
 if (make_points(domnum)) then
    call cv_i4tim_asc_lp(laps_valtime,a24time,istatus) 
    call make_fnam_lp(laps_reftime,a9time,istatus)
    call lfm_initpts(lfmprd_dir,a24time,a9time)
 endif
endif

! Write data to desired output formats.

if (out_cdf) call output_laps(lfmprd_dir,laps_data_root,domnum,laps_reftime,laps_valtime)

if (out_grib) then
   write(domnum_str,'(i2.2)') domnum
   call make_igds(proj,igds)
   call make_fnam_lp(laps_reftime,a9time,istatus)
   call make_fcst_time(laps_valtime,laps_reftime,fcst_hhmm,istatus)
   if (fcst_hhmm(3:4) == "00") fcst_hhmm(3:4)="  "
   gribfile=trim(lfmprd_dir)//'/d'//domnum_str//'/grib/'//a9time//'00'//trim(fcst_hhmm)//'.grib'
   print*,' '
   print*,'Writing 2d fields to grib: ',trim(gribfile)
   call open_grib_c(gribfile,funit)
   if (verbose) print*,'Opened funit =',funit,' for ',trim(gribfile)
   startbyte=1
   nbytes=0
   call grib_sfc_vars(laps_reftime,laps_valtime)
   print*,' '
   print*,'Writing 3d fields to grib: ',trim(gribfile)
   call grib_ua_vars(laps_reftime,laps_valtime)
   call close_grib_c(funit)
   if (verbose) print*,'Number of bytes written in grib:',nbytes
endif

call dealloc_grid('horiz')
call dealloc_grid('surface')
call dealloc_grid('isobaric')

stop ' '
900 continue
call usage

end

!===============================================================================

subroutine usage

implicit none

print*,'Usage: lfmpost.exe ''model type'' ''filename'' ''grid no'' ''laps i4time'' ''fcst time (sec)'' ''laps_data_root'''

stop

return
end
