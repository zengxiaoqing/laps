subroutine output_laps_rg(laps_data_root,mtype,domnum_in,laps_reftime,laps_valtime &
                         ,pgrid,sgrid,name2d,name3d,com2d,com3d,lvls3d,nvar2dout,nvar3dout,pres_1d,lx,ly,lz)

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

!use lfmgrid

implicit none
   
integer :: domnum_in,laps_reftime,laps_valtime,istatus,fnlen,extended
character(len=*) :: laps_data_root
character(len=256) :: output_file,donefile,command
character(len=255) :: output_dir,cdl_dir
character(len=9) :: gtime
character(len=5) :: fcst_hhmm
character(len=2) :: domnum_str
character(len=30):: mtype       

real pres_1d(lz) ! pa

!beka
integer ISTAT, I4_elapsed, ishow_timer, init_timer

integer lx,ly,lz,nvar2dout,nvar3dout

! 3D declarations
real pgrid(lx,ly,nvar3dout*lz)
real lprs(lz)
real cdl_levels(lz)
integer lvls3d(nvar3dout*lz) ! mb for pressure grid  
character*132 com3d(nvar3dout*lz)  
character*10 units3d(nvar3dout*lz)  
character*4 lvltype3d(nvar3dout*lz)  
character*3 name3d(nvar3dout*lz)    

! 2D declarations
real sgrid(lx,ly,nvar2dout)
integer lvls2d(nvar2dout)
character*132 com2d(nvar2dout)
character*10 units2d(nvar2dout)
character*4 lvltype2d(nvar2dout)
character*3 name2d(nvar2dout)

logical verbose /.true./

! LW assign domnum_in, from calling args, to domnum, declared in lfmgrid module
if (verbose) then
   print*,' '
   print*,'subroutine output_laps_rg - outputting LAPS format (fua/fsf) files...'
endif

write(6,*)' lx,ly,lz,nvar2dout,nvar3dout ',lx,ly,lz,nvar2dout,nvar3dout

cdl_dir=trim(laps_data_root)//'/cdl/'

! Lets make the fua file first (contains 3d variables).  
if (lz .gt. 1) then                  

  write(6,*)' pgrid is ',pgrid(1,1,:) 

  cdl_levels=pres_1d(lz:1:-1)*0.01

! Write out the 3D stuff using LAPS library routine

  output_dir=trim(laps_data_root)//'/lapsprd/fua/'//trim(mtype)//'/'

! Build the output file name so we can create a "donefile" if 
! running in realtime mode.

  call make_fnam_lp(laps_reftime,gtime,istatus)

  call make_fcst_time(laps_valtime,laps_reftime,fcst_hhmm,istatus)

  call cvt_fname_v3(output_dir,gtime,fcst_hhmm,'fua',3,output_file  &
                   ,fnlen,istatus)

  I4_elapsed = ishow_timer()
	
  print*,' '
  print*,'Writing 3d fields to netcdf: ',lx,ly,lz*nvar3dout,trim(output_file)        
  call system('rm -f '//trim(output_file))                 
  write(6,*)' pgrid(1,1,:) = ',pgrid(1,1,:)
  write(6,*)' lvls3d = ',lvls3d
  write(6,*)' name3d = ',name3d
  write(6,*)' com3d = ',com3d
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

endif ! Do 3d variables

! Do 2d variables.

lvls2d(:)=0.

output_dir=trim(laps_data_root)//'/lapsprd/fsf/'//trim(mtype)//'/'

call make_fnam_lp(laps_reftime,gtime,istatus)

call make_fcst_time(laps_valtime,laps_reftime,fcst_hhmm,istatus)

call cvt_fname_v3(output_dir,gtime,fcst_hhmm,'fsf',3,output_file  &
                 ,fnlen,istatus)

print*,'Writing 2d fields to netcdf: ',trim(output_file)
call system('rm -f '//trim(output_file))                 
print*,'lx,ly,nvar2dout = ',lx,ly,nvar2dout
print*,'output_dir = ',output_dir
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

return
end
