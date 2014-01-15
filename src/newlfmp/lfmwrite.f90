subroutine output_laps(lfmprd_dir,laps_data_root,domnum_in,laps_reftime,laps_valtime)

! Creates the LAPS *.fua and *.fsf files for model output.  

! Use of this routine requires a valid MM5_DATA_ROOT that contains
! a subdirectory in MM5_DATA_ROOT/mm5prd/dxx/fua (and fsf), where the
! xx is the two digit domain number.  The fsf.cdl and fua.cdl files
! in LAPS_DATA_ROOT/cdl should also have dimensions equal to lx/ly.
 
! History
! =======
! Initial version:  Brent Shaw, NOAA/FSL 5 Jan 01
! Modified to output files to domain specific directories in
! MM5_DATA_ROOT instead of LAPS_DATA_ROOT.  2 Nov 01

! Note that the only variable that is supported in the netCDF files but
! not produced by this routine is the fire index.  

use lfmgrid

implicit none
   
integer :: domnum_in,laps_reftime,laps_valtime,istatus,fnlen,extended
real, allocatable, dimension(:) :: cdl_levels
character(len=*) :: lfmprd_dir,laps_data_root
character(len=256) :: output_file,donefile
character(len=255) :: output_dir,cdl_dir
character(len=9) :: gtime
character(len=5) :: fcst_hhmm
character(len=2) :: domnum_str
logical :: outdir_defined

!beka
integer ISTAT, I4_elapsed, ishow_timer, init_timer


! LW assign domnum_in, from calling args, to domnum, declared in lfmgrid module
domnum = domnum_in

if (verbose) then
   print*,' '
   print*,'Outputting LAPS format (fua/fsf) files...'
endif

write(domnum_str,'(i2.2)') domnum
allocate(cdl_levels(lz))
cdl_levels=lprs(lz:1:-1)*0.01

! Lets make the fua file first (contains 3d variables).  
if (trim(mtype) /= 'st4') then

! Write out the 3D stuff using LAPS library routine

  write(6,*)' write_to_lapsdir = ',write_to_lapsdir

  if (.not. write_to_lapsdir) then
     output_dir=trim(lfmprd_dir)//'/d'//domnum_str//'/fua/'
  else
     output_dir=trim(laps_data_root)//'/lapsprd/fua/'//trim(mtype)//'/'
  endif
  cdl_dir=trim(laps_data_root)//'/cdl/'

! Build the output file name so we can create a "donefile" if 
! running in realtime mode.

  call make_fnam_lp(laps_reftime,gtime,istatus)

  call make_fcst_time(laps_valtime,laps_reftime,fcst_hhmm,istatus)

  call cvt_fname_v3(output_dir,gtime,fcst_hhmm,'fua',3,output_file  &
                   ,fnlen,istatus)
!beka

	I4_elapsed = ishow_timer()
	
  print*,' '
  print*,'Writing 3d fields to netcdf: ',trim(output_file)
  inquire(FILE=trim(output_dir),EXIST=outdir_defined)
  if(outdir_defined .eqv. .false.)then
      write(6,*)' ERROR: output directory does not exist'
      stop
  endif
  lvls3d=lvls3d*0.01
  call write_laps_lfm(laps_reftime,laps_valtime,trim(output_dir),trim(cdl_dir)   &
                     ,'fua'                                          &
                     ,lx,ly,lz*nvar3dout,lz*nvar3dout,name3d,lvls3d&
                     ,lvltype3d,units3d,com3d,lz,cdl_levels          &
                     ,pgrid,istatus)

  if (istatus /= 1) then
     print*,'Error writing LAPS 3D (fua) netcdf file.'
  elseif (verbose) then
     print*,'Done writing 3d netcdf file.'
  endif

!beka
	I4_elapsed = ishow_timer()

  if (realtime .and. istatus == 1 .and. make_donefile) then
     donefile=trim(output_file)//'.done'
     open(unit=2,file=donefile,status='unknown')
     close(2)
  endif
  deallocate(cdl_levels)

endif

! Do 2d variables.

lvls2d(:)=0.

if (.not. write_to_lapsdir) then
   output_dir=trim(lfmprd_dir)//'/d'//domnum_str//'/fsf/' 
else
   output_dir=trim(laps_data_root)//'/lapsprd/fsf/'//trim(mtype)//'/'
endif

cdl_dir=trim(laps_data_root)//'/cdl/'
call make_fnam_lp(laps_reftime,gtime,istatus)

call make_fcst_time(laps_valtime,laps_reftime,fcst_hhmm,istatus)

call cvt_fname_v3(output_dir,gtime,fcst_hhmm,'fsf',3,output_file  &
                 ,fnlen,istatus)

print*,'Writing 2d fields to netcdf: ',trim(output_file)
inquire(FILE=trim(output_dir),EXIST=outdir_defined)
if(outdir_defined .eqv. .false.)then
    write(6,*)' ERROR: output directory does not exist'
    stop
endif
    
call write_laps_lfm(laps_reftime,laps_valtime,trim(output_dir),cdl_dir  &
                   ,'fsf'                                               &
                   ,lx,ly,nvar2dout,nvar2dout,name2d,lvls2d,lvltype2d   &
                   ,units2d,com2d,1,0.                                  &  
                   ,sgrid,istatus)

if (istatus /= 1) then
   print*,'Error writing LAPS 2d (fsf) netcdf file.'
elseif (verbose) then
   print*,'Done writing 2d netcdf file.'
endif

if (realtime .and. istatus == 1 .and. make_donefile) then
   donefile=trim(output_file)//'.done'
   open(unit=2,file=donefile,status='unknown')
   close(2) 
endif

return
end

!===============================================================================

subroutine grib_sfc_vars(laps_reftime,laps_valtime)

use lfmgrid
use grib

implicit none

integer :: laps_reftime,laps_valtime,istatus    &
          ,id(27),yyyyr,mmr,ddr,hhr,minr,itype  &
          ,timeunit,timeperiod1,timeperiod2     &
          ,fcsttime_now,fcsttime_prev,itot      &
          ,n,shape(1)

real, pointer, dimension(:,:) :: fld2d
real, allocatable, dimension(:) :: fld

character(len=24) :: atime 
character(len=3) :: amonth,amonths(12)

data amonths/'JAN','FEB','MAR','APR','MAY','JUN'  &
            ,'JUL','AUG','SEP','OCT','NOV','DEC'/

! Compute year, month, day of month, hour, and minute from laps_reftime.

call cv_i4tim_asc_lp(laps_reftime,atime,istatus) 
read(atime,'(i2.2,x,a3,x,i4.4,x,i2.2,x,i2.2)') ddr,amonth,yyyyr,hhr,minr

do mmr=1,12
   if (amonth == amonths(mmr)) exit
enddo

! Determine appropriate timeunit.

if (mod(precip_dt,3600) == 0) then
! Time unit in hours.
   timeunit=1
   fcsttime_now=(laps_valtime-laps_reftime)/3600
   if (fcsttime_now > 0) then
      fcsttime_prev=fcsttime_now-(precip_dt/3600)
   else
      fcsttime_prev=0
   endif
else
! Time unit in minutes.
   timeunit=0
   fcsttime_now=(laps_valtime-laps_reftime)/60
   if (fcsttime_now > 0) then
      fcsttime_prev=fcsttime_now-(precip_dt/60)
   else
      fcsttime_prev=0
   endif
endif
itype=0
    
! Grib up each variable.

allocate(fld(lx*ly))
shape(1)=lx*ly
do n=1,nvar2dout
   if (gribit(n)) then
      fld2d=>sgrid(1:lx,1:ly,n)
      if (minval(fld2d) >= rmsg) cycle
      fld=reshape(fld2d,shape)
      if (verbose) write(6,'('' Gribbing '',a25,'', Min/Max = '',2(f10.4,1x))')  &
                        com2d(n)(1:25),minval(fld),maxval(fld)
      if (timerange(n) == 0) then 
         timeperiod1=fcsttime_now
         timeperiod2=0
      elseif (timerange(n) == 4) then
         timeperiod1=fcsttime_prev
         timeperiod2=fcsttime_now
      endif
      call make_id(table_version,center_id,subcenter_id,process_id          &
                  ,param(n),leveltype(n),level1(n),level2(n),yyyyr,mmr,ddr  &
                  ,hhr,minr,timeunit,timerange(n),timeperiod1,timeperiod2   &
                  ,scalep10(n),id) 
      call write_grib(itype,fld,id,igds,funit,startbyte,itot,istatus)
      nbytes=nbytes+itot
      startbyte=nbytes+1 
   endif
enddo
deallocate(fld)

return
end

!===============================================================================

subroutine grib_ua_vars(laps_reftime,laps_valtime)

use lfmgrid
use grib

implicit none

integer :: laps_reftime,laps_valtime,istatus  &
          ,id(27),yyyyr,mmr,ddr,hhr,minr      &
          ,itype,lvltype,lvl1,lvl2,tmrange     &
          ,timeunit,timeperiod1,timeperiod2   &
          ,fcsttime_now,fcsttime_prev,itot    &
          ,k,n,shape(1)

real, pointer, dimension(:,:) :: fld2d
real, allocatable, dimension(:) :: fld

character(len=24) :: atime 
character(len=3) :: amonth,amonths(12)

data amonths/'JAN','FEB','MAR','APR','MAY','JUN'  &
            ,'JUL','AUG','SEP','OCT','NOV','DEC'/

! Compute year, month, day of month, hour, and minute from laps_reftime

call cv_i4tim_asc_lp(laps_reftime,atime,istatus)
read(atime,'(i2.2,x,a3,x,i4.4,x,i2.2,x,i2.2)') ddr,amonth,yyyyr,hhr,minr

do mmr=1,12
   if (amonth == amonths(mmr)) exit
enddo

! Determine appropriate timeunit.

if (mod(precip_dt,3600) == 0) then
! Time unit in hours.
   timeunit=1
   fcsttime_now=(laps_valtime-laps_reftime)/3600
   if (fcsttime_now > 0) then
      fcsttime_prev=fcsttime_now-(precip_dt/3600)
   else
      fcsttime_prev=0
   endif
else
! Time unit in minutes.
   timeunit=0
   fcsttime_now=(laps_valtime-laps_reftime)/60
   if (fcsttime_now > 0) then
      fcsttime_prev=fcsttime_now-(precip_dt/60)
   else
      fcsttime_prev=0
   endif
endif

! Grib up each variable at each level.

itype=0
lvltype=100
lvl2=0
tmrange=0
timeperiod1=fcsttime_now
timeperiod2=0

allocate(fld(lx*ly))
shape(1)=lx*ly
do k=1,lz
   lvl1=nint(lprs(k))*0.01
   if (lvl1 <= 1000) then
      do n=1,nvar3dout
         if (gribitua(n)) then
            fld2d=>pgrid(1:lx,1:ly,(n-1)*lz+k)
            if (minval(fld2d) >= rmsg) cycle
            fld=reshape(fld2d,shape)
            if (verbose) write(6,'('' Gribbing '',a25,'', level='',i4,''mb, Min/Max = '',2(f10.4,1x))')  &
                              com3d((n-1)*lz+k)(1:25),lvl1,minval(fld),maxval(fld)
            call make_id(table_version,center_id,subcenter_id,process_id      &
                        ,paramua(n),lvltype,lvl1,lvl2,yyyyr,mmr,ddr     &
                        ,hhr,minr,timeunit,tmrange,timeperiod1,timeperiod2  &
                        ,scalep10(n),id) 
            call write_grib(itype,fld,id,igds,funit,startbyte,itot,istatus)
            nbytes=nbytes+itot
            startbyte=nbytes+1 
         endif
      enddo
   endif
enddo
deallocate(fld)

return
end
