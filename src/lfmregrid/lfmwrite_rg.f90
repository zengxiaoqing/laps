subroutine output_laps_rg(lfmprd_dir,laps_data_root,domnum_in,laps_reftime,laps_valtime &
                         ,sgrid,name2d,com2d,nvar2dout,nvar3dout,lx,ly,lz)

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
character(len=*) :: lfmprd_dir,laps_data_root
character(len=256) :: output_file,donefile
character(len=255) :: output_dir,cdl_dir
character(len=9) :: gtime
character(len=5) :: fcst_hhmm
character(len=2) :: domnum_str

!beka
integer ISTAT, I4_elapsed, ishow_timer, init_timer

integer lx,ly,lz,nvar2dout,nvar3dout

! 3D declarations
real pgrid(lx,ly,nvar3dout*lz)
real lprs(lz)
real cdl_levels(lz)
integer lvls3d(lz)   
character*132 com3d(lz)  
character*10 units3d(lz)  
character*4 lvltype3d(lz)  
character*3 name3d(lz)    

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
   print*,'Outputting LAPS format (fua/fsf) files...'
endif

cdl_levels=lprs(lz:1:-1)*0.01

! Lets make the fua file first (contains 3d variables).  
if (lz .gt. 1) then                  

! Write out the 3D stuff using LAPS library routine

! output_dir=trim(laps_data_root)//'/lapsprd/fua/'//trim(mtype)//'/'
  output_dir=trim(lfmprd_dir)//'/'                                                 
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

endif ! Do 3d variables

! Do 2d variables.

lvls2d(:)=0.

!output_dir=trim(laps_data_root)//'/lapsprd/fsf/'//trim(mtype)//'/'
output_dir=trim(lfmprd_dir)//'/'                                                 

cdl_dir=trim(laps_data_root)//'/cdl/'
call make_fnam_lp(laps_reftime,gtime,istatus)

call make_fcst_time(laps_valtime,laps_reftime,fcst_hhmm,istatus)

call cvt_fname_v3(output_dir,gtime,fcst_hhmm,'fsf',3,output_file  &
                 ,fnlen,istatus)

print*,'Writing 2d fields to netcdf: ',trim(output_file)
    
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
