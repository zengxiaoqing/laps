!====================================================================
!  This wrapper is designed to scan all possible airborne radar data
!  files and to read all superobs data.
!
!  History:
!
!       Creation: Yuanfu Xie August 2009.
!====================================================================

subroutine wrapper(i4time,halftime,path,maxdat,numdat,obsi4t, &
                   lat,lon,hgh,azi,elv,dis,wnd)

  implicit none

  character*256, intent(in) :: path	! Path to the airborne data
  integer, intent(in) :: i4time,halftime,maxdat
  integer, intent(out) ::numdat
  integer, intent(out) :: obsi4t(maxdat)	! obs i4time
  real,    intent(out) :: lat(maxdat),lon(maxdat),hgh(maxdat), &
                          azi(maxdat),elv(maxdat) ! Azimuth/elevation angles
  real,    intent(out) :: dis(maxdat),wnd(maxdat) ! Radius and radial wind

  ! Local variables:
  integer, parameter :: max_files = 2000
  integer :: ifile,nfiles,l1,l2,istatus
  integer :: beginI4,endngI4
  character*256 :: c_filenames(max_files),filter
  character*13 :: beginTime,endngTime,obsTime
  character*9 :: begin9,endng9,wfo_fname13_to_fname9

  ! Initialize numdat:
  numdat = 0

  ! List of all available files:
  ! call get_file_time(path,nfiles,c_filenames,max_files,istatus)
  filter = '*'
  call getfilenames_c(path,c_filenames,nfiles,filter,istatus)

  ! For file in the time window:
  ! Hardcopy for handling time after year 2000 only: !!!!
  beginTime(1:2) = '20'
  endngTime(1:2) = '20'
  beginTime(9:9) = '_'
  endngTime(9:9) = '_'
  do ifile=1,nfiles
    ! YYMMDD:
    beginTime(3:8) = c_filenames(ifile)(1:6)
    endngTime(3:8) = c_filenames(ifile)(1:6)
    ! HHMM:
    beginTime(10:13) = c_filenames(ifile)(10:13)
    endngTime(10:13) = c_filenames(ifile)(15:18)

    ! Convert to fname9:
    begin9 = wfo_fname13_to_fname9(beginTime)
    endng9 = wfo_fname13_to_fname9(endngTime)

    ! Convert to i4 time:
    call cv_asc_i4time(begin9,beginI4,istatus)
    call cv_asc_i4time(endng9,endngI4,istatus)

    ! Find a overlap:
    if (i4time+halftime .ge. beginI4 .and. i4time-halftime .le. endngI4) then

      ! Open data file:
      call s_len(path,l1)
      call s_len(c_filenames(ifile),l2)
      open(13,file=path(1:l1)//c_filenames(ifile)(1:l2),status='old')

      ! Read a record:
      numdat = numdat+1
      if (numdat .gt. maxdat) then
        print*,'Error: too many radial wind data: ',numdat,' > ',maxdat
        stop
      endif

 1    read(13,*,end=2) obsTime(1:12),lat(numdat),lon(numdat),hgh(numdat), &
                       azi(numdat),elv(numdat),dis(numdat),wnd(numdat)

      obsTime(10:13) = obsTime(9:12)
      obsTime(9:9) = '_'
      obsTime(1:9) = wfo_fname13_to_fname9(obsTime(1:13))

      ! obs i4time:
      call cv_asc_i4time(obsTime(1:9),obsi4t(numdat),istatus)
      ! Skip invalid obs time:
      if (obsi4t(numdat) .lt. i4time-halftime .or. &
          obsi4t(numdat) .gt. i4time+halftime) then
        numdat = numdat-1
      endif
      ! Next record:
      goto 1
 2    close(13)

    endif
  enddo

end subroutine wrapper
