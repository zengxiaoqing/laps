subroutine read_laps_radar
!==================================================================
!  This routine reads in radar through LAPS ingest codes.
!
!  Note:
!  1. The original RDLAPSRDR uses LAPS get_multiradar_vel to ingest 
!     which requires the number of gridpoints times the number of
!     radar and causes problems for large domain data. The routine
!     avoids those arrays for efficiency in computer memory by 
!     repeating what getradar.f under lib does.
!  2. This routine assumes to be used as a member of module of the
!     readobserves in readobserves.f90 and it uses all variables
!     defined in STMAS.
!
!  History: August 2012
!  Author : Yuanfu Xie at ESRL/GSD/FAB
!==================================================================

  use prmtrs_stmas

  implicit none

  ! Get radar radial wind obs:
  ! print*,'Reading radial wind... DUE TO RADAR FORMAT CHANGE: NO RADIAL WIND ingested!' 
  if (fcstgrd(1) .lt. 10000) call radialwind

  ! Get radar reflectivity obs:
  print*,'Reading reflectivity...'
  call reflectivity

end subroutine read_laps_radar

subroutine radialwind
!==================================================================
!  This routine reads in radial wind data through LAPS ingest.
!
!  History: August 2012
!  Author : Yuanfu Xie at ESRL/GSD/FAB
!==================================================================

  use prmtrs_stmas
  use read_backgrd
  use readobserves
  use mem_namelist

  implicit none

  include 'main_sub.inc'

  ! Local variables:
  character :: ext*31, rname*4, sid*8, dir*300
  integer   :: i,j,k,ll,iradar,len_dir,istatus,istat,i4radar,num
  logical   :: apply_map ! LAPS does not use this for radial wind
  real      :: volnyqt,rlat,rlon,rhgt,op(4),az,sr,ea

  ! For radial wind data on maxgrid:
  INTEGER   :: istart,iend,jstart,jend
  REAL      :: ri,rj,grid_spacing_actual_mx,grid_spacing_actual_my
  REAL      :: xradius,yradius

  ! Functions:
  integer :: init_timer,ishow_timer
  real :: height_of_level
  
  ! Allocatable arrays:
  real,allocatable :: velo(:,:,:), nyqt(:,:,:),zero(:,:,:)

  ! Allocate memory of local arrays:
  allocate(velo(fcstgrd(1),fcstgrd(2),fcstgrd(3)), &
           nyqt(fcstgrd(1),fcstgrd(2),fcstgrd(3)), &
           zero(fcstgrd(1),fcstgrd(2),fcstgrd(3)), &
           stat=istatus)

  ! Ingest:
  PRINT*,'Starting ingest radial wind: ', init_timer()
  do iradar=1,max_radars ! look through all radar

    write(6,*)' Looking for potential radar # ',iradar

    ! Get file:
    if (iradar .lt. 10)                         write(ext,1) iradar
    if (iradar .ge. 10  .and. iradar .lt.  100) write(ext,2) iradar
    if (iradar .ge. 100 .and. iradar .lt. 1000) write(ext,3) iradar
    if (iradar .ge. 1000) then
      print*,'Too many radars: modify this code and rerun!'
    endif
1   format('v0',i1,' ')
2   format('v', i2,' ')
3   format('v', i3,' ')

    ! Get that directory:
    call get_directory(ext,dir,len_dir)
    dir = dir(1:len_dir)//'*.'//ext

    ! Frequency of reading radial wind data:
    do ll=2,2 ! read data at current ONLY !and previous cycle

      call get_file_time(dir,lapsi4t+(ll-2)*INT(0.5*icycle),i4radar)

      ! Check file availability:
      if (abs(i4radar-(lapsi4t+(ll-2)*0.5*icycle)) .ge. 0.5*icycle) then
        write(6,10) lapsi4t+(ll-2)*icycle,dir(1:len_dir)
10      format('No radial wind avail at ',i10,' under ', &
          /,' in directory: ',a50)
      else ! find good file:
        num = 0 ! number of grid points where radial wind has value
        call read_radar_vel(i4radar,apply_map, &
          fcstgrd(1),fcstgrd(2),fcstgrd(3),ext, &
          velo,nyqt,volnyqt,rlat,rlon,rhgt,rname,num,istatus,istat)
        if (num .gt. 0) then
          ! radar QC for unfolding radial wind:
          zero = 0.0
          call qc_radar_obs(fcstgrd(1),fcstgrd(2),fcstgrd(3), &
                            rmissing,fcstgrd(1),fcstgrd(2),0,0, &
                            velo,nyqt,num,latitude,longitud, &
                            rlat,rlon,rhgt,zero,zero, &
                            bk0(1,1,1,ll,1),bk0(1,1,1,ll,2),volnyqt, &
                            l_correct_unfolding,l_grid_north,istatus)

          ! Estimate this radar coverage:
          CALL latlon_to_rlapsgrid(rlat,rlon,latitude,longitud, &
                                   fcstgrd(1),fcstgrd(2),ri,rj,istatus)
 
          ! Grid spacing:
          CALL get_grid_spacing_actual_xy(rlat,rlon &
                            ,grid_spacing_actual_mx &
                            ,grid_spacing_actual_my &
                            ,istatus)

          ! Radius on fcstgrid:
          xradius = 300000.0/grid_spacing_actual_mx  ! Radius of radial wind in X
          yradius = 300000.0/grid_spacing_actual_my  ! Radius of radial wind in Y

          ! Start and end points on fcstgrd:
          istart = MAX(NINT(ri-xradius),1)
          jstart = MAX(NINT(rj-yradius),1)
          iend = MIN(NINT(ri+xradius),fcstgrd(1))
          jend = MIN(NINT(rj+yradius),fcstgrd(2))

          ! Pass radial wind to STMAS data structure:
          do k=1,fcstgrd(3)
          do j=jstart,jend
          do i=istart,iend

            if (velo(i,j,k) .ne. rmissing) then

              ! Calculate azimuth/elevation/slant range:
              call latlon_to_radar(latitude(i,j),longitud(i,j), &
                                   bk0(i,j,k,ll,3),az,sr,ea,rlat,rlon,rhgt)

              op(1) = longitud(i,j)
              op(2) = latitude(i,j)
              op(3) = z_fcstgd(k)
              ! Treat data as observed at analysis time frames:
              op(4) = 0.5*icycle ! i4radar-itime2(1)

              call handleobs_sigma(op,velo(i,j,k),0.5, & ! radial obs error: 0.5
                                   numstat+1,nallobs,1,az,ea,sid)
            endif
          enddo
          enddo
          enddo
        endif
      endif

    enddo

  enddo
  PRINT*,'Ending ingest radial wind: ', ishow_timer()

  ! Release local allocatables:
  deallocate(velo,nyqt,zero,stat=istatus)

end subroutine radialwind

subroutine reflectivity
!==================================================================
!  This routine reads in reflectivity data through LAPS ingest for
!  deriving humidity (SH) lower bounds.
!  Note: cloud fraction is also read in for detecting non-cloudness.
!
!  History: August 2012
!  Author : Yuanfu Xie at ESRL/GSD/FAB
!==================================================================

  use prmtrs_stmas
  use read_backgrd
  use readobserves
  use mem_namelist

  implicit none

  ! Local variables:
  character :: ext*31, rname*4, unit*10, comment*150
  integer :: i,j,k,ll,iqc,nref,n2d,n3d,istatus,ix0,ix1,iy0,iy1
  integer :: istatus2d(fcstgrd(1),fcstgrd(2)), &
             istatus3d(fcstgrd(1),fcstgrd(2))
  logical :: reflectivity_bound
  real    :: rhc,tref ! LAPS uses these scalars
  real    :: rmax(2),rlow(2)

  ! Functions:
  real :: make_ssh

  ! Include LAPS cloud parameters:
  include 'laps_cloud.inc'

  ! Allocatables: Reflectivity, cloud fraction on p-coor and on height-coord:
  real, allocatable :: refl(:,:,:,:),cldf(:,:,:,:),cldh(:,:,:,:)
  real :: cld_hgt(kcloud),cld_prs(kcloud)

  ! Allocate memory:
  allocate(refl(fcstgrd(1),fcstgrd(2),fcstgrd(3),3), &
           cldf(fcstgrd(1),fcstgrd(2),fcstgrd(3),3), &
           cldh(fcstgrd(1),fcstgrd(2),kcloud,3), &
    stat=istatus)
  
  ! Reference temperaure for converison between SH and RH:
  tref = 0.0 !-132.0

  ! Get reflectivity at pressure levels:
  ext = 'lps '
  refl = 0.0
  do ll=1,2
    call get_laps_3dgrid(lapsi4t+(ll-2)*icycle,icycle/2,i, &
                         fcstgrd(1),fcstgrd(2),fcstgrd(3), &
                         ext,'ref',unit,comment,refl(1,1,1,ll),istatus)
    if (i .ne. lapsi4t+(ll-2)*icycle) then
      print*,'Does not find a lps at: ',lapsi4t+(ll-2)*icycle
      cycle
    endif
  enddo
  ! Interpolate: 
  refl(1:fcstgrd(1),1:fcstgrd(2),1:fcstgrd(3),1) = &
    0.5*(refl(1:fcstgrd(1),1:fcstgrd(2),1:fcstgrd(3),1)+ &
         refl(1:fcstgrd(1),1:fcstgrd(2),1:fcstgrd(3),2))
  refl(1:fcstgrd(1),1:fcstgrd(2),1:fcstgrd(3),3) = &
    2.0*refl(1:fcstgrd(1),1:fcstgrd(2),1:fcstgrd(3),2)- &
        refl(1:fcstgrd(1),1:fcstgrd(2),1:fcstgrd(3),1) ! extrapolation
  
  ! Reflectivity derived bounds:
  bk0(1:fcstgrd(1),1:fcstgrd(2),1:fcstgrd(3),1:fcstgrd(4),numstat+1) = 0.0
  rmax(1) = 45.0
  rlow(1) = 0.5
  rmax(2) = 30.0
  rlow(2) = 0.6
  do ll=1,fcstgrd(4)
  do k=1,fcstgrd(3)
  do j=1,fcstgrd(2)
  do i=1,fcstgrd(1)
    if (refl(i,j,k,ll) .gt. 5.0 .and. refl(i,j,k,ll) .le. 100) then

      ! iqc: 1 rain; 2 snow:
      iqc = 2
      if (bk0(i,j,k,ll,temprtur) .gt. 273.15) iqc = 1

      if (refl(i,j,k,ll) .gt. rmax(iqc)) then
        bk0(i,j,k,ll,numstat+1) = &
          make_ssh(z_fcstgd(k)/100.0,bk0(i,j,k,ll,temprtur)-273.15,1.0,tref)
      else
        bk0(i,j,k,ll,numstat+1) = &
          make_ssh(z_fcstgd(k)/100.0,bk0(i,j,k,ll,temprtur)-273.15, &
            1.0,tref)
  !          rlow(iqc)+(1.0-rlow(iqc))*refl(i,j,k,ll)/rmax(iqc),tref)
      endif
    endif
  enddo
  enddo
  enddo
  enddo

  ! Read LAPS cloud fraction on LAPS cld_hgt coordinate and then interpolate to pressure levels:
  ext = 'lc3 '
  cldf = 0.0
  do ll=1,2 ! lc3 at current and previous time frames
    call get_clouds_3dgrid(lapsi4t+(ll-2)*icycle,i, &
                         fcstgrd(1),fcstgrd(2),kcloud, &
                         ext,cldh(1,1,1,ll),cld_hgt,cld_prs,istatus)
    if (i .ne. lapsi4t+(ll-2)*icycle) then
      print*,'Does not find a lc3 at: ',lapsi4t+(ll-2)*icycle
      cycle
    endif

    ! Interpolate from LAPS_height to analysis pressure coordinate:
    print*,'Interpolating lc3 to lcp...'
    call interp_height_pres_fast(fcstgrd(1),fcstgrd(2),fcstgrd(3),kcloud,cldf(1,1,1,ll), &
         cldh(1,1,1,ll),bk0(1,1,1,ll,pressure),cld_hgt,istatus)

  enddo
  print*,'Cloud fraction before interpolation: ',maxval(cldf(:,:,:,1)),minval(cldf(:,:,:,1))

  ! Interpolate cloud fraction to all time frames:
  cldf(1:fcstgrd(1),1:fcstgrd(2),1:fcstgrd(3),2) = &
    0.5*(cldf(1:fcstgrd(1),1:fcstgrd(2),1:fcstgrd(3),1)+ &
         cldf(1:fcstgrd(1),1:fcstgrd(2),1:fcstgrd(3),2))
  cldf(1:fcstgrd(1),1:fcstgrd(2),1:fcstgrd(3),3) = &
    2.0*cldf(1:fcstgrd(1),1:fcstgrd(2),1:fcstgrd(3),2)- &
        cldf(1:fcstgrd(1),1:fcstgrd(2),1:fcstgrd(3),1) ! extrapolation
  
  ! Adjust bound based on cloud:
  reflectivity_bound = .false.
  do ll=1,fcstgrd(4)
  do k=1,fcstgrd(3)
  do j=1,fcstgrd(2)
  do i=1,fcstgrd(1)
    ! Temporarily remove sh lower bounds where there is no cloud:
    if (cldf(i,j,k,ll) .lt. 0.1 .and. refl(i,j,k,ll) .gt. 5.0) &
      bk0(i,j,k,ll,numstat+1) = 0.0

    ! Covering the following 5 lines: Testing no bound for cloudy area without reflectivity:
    !if (cldf(i,j,k,ll) .ge. 0.1) then
    !  rhc = make_ssh(z_fcstgd(k)/100.0,bk0(i,j,k,ll,temprtur)-273.15, &
    !                 cldf(i,j,k,ll)**0.2,0.0)
    !  bk0(i,j,k,ll,numstat+1) = amax1(bk0(i,j,k,ll,numstat+1),rhc)
    !endif

    ! RH =100% if both cloud and reflectivity occur:
    if (cldf(i,j,k,ll) .ge. 0.1 .and. refl(i,j,k,ll) .ge. 5.0) &
      bk0(i,j,k,ll,numstat+1) = &
          make_ssh(z_fcstgd(k)/100.0,bk0(i,j,k,ll,temprtur)-273.15,1.0,tref)
    if (ll .eq. 2 .and. cldf(i,j,k,ll) .ge. 0.1 .and. refl(i,j,k,ll) .ge. 5.0) &
      reflectivity_bound = .true.

  enddo
  enddo
  enddo
  enddo
  print*,'Reflectivity bound derived: ',reflectivity_bound

  ! Assign reflectivity derived bounds to GRDBKGD0:
  do j=1,maxgrid(2)
    iy0 = float(j-1)/float(maxgrid(2)-1)*(fcstgrd(2)-1)+1
    do i=1,maxgrid(1)
      ix0 = float(i-1)/float(maxgrid(1)-1)*(fcstgrd(1)-1)+1
      do ll=1,maxgrid(4) 
      do k=1,maxgrid(3) 
        grdbkgd0(i,j,k,ll,numstat+1) = bk0(ix0,iy0,k,ll,numstat+1)
        grdbkgd0(i,j,k,ll,numstat+2) = 1000.0 ! no uppper bound now
      enddo
      enddo
  enddo
  enddo
  print*,'Max dBZ over finest grid: ',maxval(grdbkgd0(:,:,:,:,numstat+1)), &
                                      minval(grdbkgd0(:,:,:,:,numstat+1))
  ! Deallocate:
  
  deallocate(refl,cldf,cldh,stat=istatus)

end subroutine reflectivity

