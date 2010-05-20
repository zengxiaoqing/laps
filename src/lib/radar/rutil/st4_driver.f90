  !==============================================================================================!
  ! Subroutine st4_driver(infile,time,debug,err) by Brad Beechler 2010/03                        !
  !==============================================================================================!
  !                                                                                              !
  !      Purpose  - Reads in a netCDF file called <filename> that contains stage iv data         !
  !                                                                                              !
  !      Inputs   - infile  : a text string with the input file's path and name                  !
  !                 time    : the laps i4time of the analysis (needed for header info)           ! 
  !                 debug   : the level of detail to output                                      !
  !                          0 - No output (except for errors)                                   !
  !                          1 - Minimal output                                                  !
  !                          2 - Detailed output                                                 !
  !                                                                                              !
  !      Outputs  - err: a code that gives the exit status of the routine                        !
  !                     0 - Successful (passed all internal tests)                               !
  !                     1 - Couln't open the stage iv netCDF file (error from netCDF follows)    !
  !                     2 - Couldn't define LAPS domain                                          !
  !                     3 - Problem writing the output file                                      !
  !                                                                                              !
  !      Requires - libraries  :                                                                 !
  !                     netCDF                                                                   !
  !                 subroutines:                                                                 !
  !                     cdf(status, err)       error handling                                    !
  !                     alloc_err(errorcode)   error handling                                    !
  !----------------------------------------------------------------------------------------------!
  SUBROUTINE st4_driver(infile,lapsi4t,debug,err)
    USE netcdf
    IMPLICIT NONE
    Real, Parameter       :: st4_missing = 1E20                ! Missing value for stage IV
    Integer, Intent(In)   :: debug                             ! Input
    Integer*4, Intent(In) :: lapsi4t                           ! Input
    Integer, Intent(Out)  :: err                               ! Output
    Character*(*) :: infile                                    ! I/O Strings
    Integer       :: ncid, dimxid, dimyid                      ! netCDF IDs
    Integer       :: latvarid, lonvarid, datavarid             ! netCDF IDs
    Integer       :: xcount, ycount                            ! Looping
    Integer       :: num_good                                  ! Counting
    Integer       :: inxdim, inydim, outxdim, outydim, cdfout  ! Holders
    Real          :: laps_missing                              ! Missing value for LAPS
    Real          :: data_min, data_max, data_mean             ! Holders
    Real          :: la1, lo1, lov, in_dx, in_dy               ! Holders
    Real          :: lat_corners(4), lon_corners(4)            ! Holders
    Real, Allocatable, Dimension(:,:) :: indata, outdata       ! Data Arrays
    Real, Allocatable, Dimension(:,:) :: inlat, outlat         ! Data Arrays
    Real, Allocatable, Dimension(:,:) :: inlon, outlon         ! Data Arrays
    ! These are specific to LAPS routines
    Character*150 :: directory
    Character*125 :: comment
    Character*31  :: extension
    Character*10  :: units
    Character*3   :: variable
    Character*4   :: lvl_coord
    Integer       :: lvl
    Character*132 :: cmodel
    Real, Allocatable, Dimension(:,:) :: grx, gry
    Integer nxc,nyc,nzc
    Real sw(2),ne(2),rota,lat0,lon0
    Real xmin,ymin,dx,dy
    common /psgrid/nxc,nyc,nzc,lat0,lon0,rota,sw,ne
    common /pscorner/xmin,ymin,dx,dy

    lvl = 0
    err = 0
    if (debug.ne.0) write(6,'(a)') 'st4_driver v1.0 by Brad Beechler'
   
    ! This part is the netCDF side of the story where we open the stage IV
    ! file and read in all of the domain info and data
    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    CALL cdf(nf90_open(infile, nf90_nowrite, ncid), err) ! Open the file read-only.
    If (err.ne.0) RETURN
    ! Get the x and y dimensions (for some reason NCL has reversed these, no big deal)
    CALL cdf(nf90_inq_dimid(ncid, 'g5_x_0', dimyid), err) ! Get the dimid of the data
    CALL cdf(nf90_inq_dimid(ncid, 'g5_y_1', dimxid), err)
    CALL cdf(nf90_inquire_dimension(ncid, dimxid, len = inxdim), err)
    CALL cdf(nf90_inquire_dimension(ncid, dimyid, len = inydim), err)
    If (err.ne.0) RETURN ! A strange error occurred, maybe corrupt file

    Allocate(indata(inxdim,inydim), stat=err) ; If (err.ne.0) CALL alloc_error(err)
    Allocate(inlat(inxdim,inydim), stat=err) ; If (err.ne.0) CALL alloc_error(err)
    Allocate(inlon(inxdim,inydim), stat=err) ; If (err.ne.0) CALL alloc_error(err)
    
    ! Get the domain information
    CALL cdf(nf90_inq_varid(ncid, 'g5_lat_0', latvarid), err)
    CALL cdf(nf90_inq_varid(ncid, 'g5_lon_1', lonvarid), err)
    CALL cdf(nf90_get_att(ncid, latvarid, 'La1', la1), err)
    CALL cdf(nf90_get_att(ncid, latvarid, 'Lo1', lo1), err)
    CALL cdf(nf90_get_att(ncid, latvarid, 'Lov', lov), err)
    CALL cdf(nf90_get_att(ncid, latvarid, 'Dx',  in_dx), err)
    CALL cdf(nf90_get_att(ncid, latvarid, 'Dy',  in_dy), err)
    CALL cdf(nf90_get_att(ncid, latvarid, 'corners',  lat_corners), err)
    CALL cdf(nf90_get_att(ncid, lonvarid, 'corners',  lon_corners), err)
    If (err.ne.0) RETURN ! A strange error occurred, maybe corrupt file
    
    ! Get the data information this is different for different accumulation times so
    ! we will look for all three kinds, only one can be there
    cdfout = nf90_inq_varid(ncid, 'A_PCP_GDS5_SFC_acc24h', datavarid)
    If (cdfout.ne.0) cdfout = nf90_inq_varid(ncid, 'A_PCP_GDS5_SFC_acc6h',  datavarid)
    If (cdfout.ne.0) cdfout = nf90_inq_varid(ncid, 'A_PCP_GDS5_SFC_acc1h',  datavarid)
    If (cdfout.ne.0) Then ! Something is wrong and we need to abort
      err = 1 ! Set output err value to the netCDF problem code
      If (debug.ne.0) write(6,'(a)') 'Aborting from netCDF: '//nf90_strerror(cdfout)
      RETURN
    EndIf
    CALL cdf(nf90_get_var(ncid, datavarid, indata, start = (/ 1, 1 /) ), err)

    ! Get latitudes and longitudes
    CALL cdf(nf90_inq_varid(ncid, 'g5_lat_0', datavarid), err)    
    CALL cdf(nf90_get_var(ncid, datavarid, inlat, start = (/ 1, 1 /) ), err)
    CALL cdf(nf90_inq_varid(ncid, 'g5_lon_1', datavarid), err)    
    CALL cdf(nf90_get_var(ncid, datavarid, inlon, start = (/ 1, 1 /) ), err)
    
    ! Convert to meters and set missing values to LAPS values
    CALL get_r_missing_data(laps_missing,err)
    If (err.ne.1) Then
       err = 2 ! Set output err value to the LAPS problem code
       write(6,'(a)') 'Error getting the laps missing data value'
       RETURN
    EndIf
    data_min = st4_missing
    data_max = -st4_missing
    num_good = 0
    data_mean = 0
    Do xcount = 1,inxdim ! Loop over x and y
    Do ycount = 1,inydim
      If (indata(xcount,ycount).ge.st4_missing) Then
        indata(xcount,ycount) = laps_missing
      Else If (debug.gt.1) Then ! Only need to do this for debug output
        indata(xcount,ycount) = indata(xcount,ycount) / 1000.0
        num_good = num_good + 1
        If (indata(xcount,ycount).gt.data_max) data_max = indata(xcount,ycount)
        If (indata(xcount,ycount).lt.data_min) data_min = indata(xcount,ycount)
        data_mean = ((data_mean * (num_good-1)) + indata(xcount,ycount)) / num_good
      EndIf
    EndDo
    EndDo ! Loop over x and y
    
    CALL cdf(nf90_close(ncid), err) ! Close the stage IV file
    
    If (debug.gt.1) Then ! Output grid information
      write(6,'(a)') ''
      write(6,'(a)') 'Input information  '
      write(6,'(a)') '^^^^^^^^^^^^^^^^^^ '
      write(6,'(a,i4,a,i4)') 'Dimensions : ',inxdim,' x ',inydim
      !write(6,'(a,f7.2,a,f7.2,a,f7.2)') 'La1/Lo1 Lov: ',la1,' / ', lo1, '   ', lov
      !write(6,'(a,f9.2,a,f9.2)') 'Dx/Dy      : ',in_dx,' / ', in_dy
      write(6,'(a,f7.2,a,f7.2,a,f7.2,a,f7.2)') 'Lat Corners: ',lat_corners(1),', ', &
                       lat_corners(2),', ',lat_corners(3),', ',lat_corners(4)
      write(6,'(a,f7.2,a,f7.2,a,f7.2,a,f7.2)') 'Lon Corners: ',lon_corners(1),', ', &
                       lon_corners(2),', ',lon_corners(3),', ',lon_corners(4)
      write(6,'(a,i7,a,i7,a,f6.2,a)') 'Number of good values/total (%) : ', num_good, &
              ' / ',inxdim*inydim,'  (', (num_good*100.0)/(inxdim*inydim),'% )'
      write(6,'(a,f6.3,a,f6.3,a,f7.3)') 'Data in mm (w/o missing) min/mean/max : ', &
                       1000.0*data_min,' / ',1000.0*data_mean,' / ',1000.0*data_max
    EndIf

    ! This part is the laps side of things that gets the output information and writes
    ! the data to the laps grid defined in $LAPS_DATA_ROOT
    ! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    CALL get_grid_dim_xy(outxdim,outydim,err)
    If (err.ne.1) Then
      err = 2 ! Set output err value to the LAPS problem code
      If (debug.ne.0) write(6,'(a)') 'From LAPS: error reading nest7grid.parms'
      RETURN
    EndIf

    Allocate(outdata(outxdim,outydim),stat=err) ; If (err.ne.0) CALL alloc_error(err)
    Allocate(outlat(outxdim,outydim),stat=err) ; If (err.ne.0) CALL alloc_error(err)
    Allocate(outlon(outxdim,outydim),stat=err) ; If (err.ne.0) CALL alloc_error(err)
    Allocate(grx(outxdim,outydim),stat=err) ; If (err.ne.0) CALL alloc_error(err)
    Allocate(gry(outxdim,outydim),stat=err) ; If (err.ne.0) CALL alloc_error(err)
    
    CALL read_static_grid(outxdim,outydim,'LAT',outlat,err)
    If (err.ne.1) Then
      err = 2
      If (debug.ne.0) write(6,'(a)') 'From LAPS: error reading latitude'
      RETURN
    EndIf
    CALL read_static_grid(outxdim,outydim,'LON',outlon,err)
    If (err.ne.1) Then
      err = 2
      If (debug.ne.0) write(6,'(a)') 'From LAPS: error reading longitude'
      RETURN
    EndIf
    
    ! Set global values for LAPS
    nxc   = inxdim
    nyc   = inydim
    nzc   = 1
    lat0  = 90.0
    lon0  = lov
    rota  = 0.0
    sw(1) = lat_corners(1) ! Southwest
    sw(2) = lon_corners(1)
    ne(1) = lat_corners(3) ! Northeast
    ne(2) = lon_corners(3)
    xmin  = 1
    ymin  = 1
    dx    = in_dx
    dy    = in_dy

    ! Remap the data to the LAPS domain
    CALL init_hinterp(inxdim,inydim,outxdim,outydim,'PS',outlat,outlon,grx,gry,0,cmodel,.false.)
    CALL hinterp_field(inxdim,inydim,outxdim,outydim,1,grx,gry,indata,outdata,.false.)

    ! Quality control negative values produced by the cubic interpolation 
    Do xcount = 1,outxdim ! Loop over x and y 
    Do ycount = 1,outydim
      If (outdata(xcount,ycount).lt.0.0) outdata(xcount,ycount) = 0.0
    EndDo
    EndDo ! Loop over x and y

    data_min = laps_missing
    data_max = -laps_missing
    num_good = 0
    data_mean = 0
    Do xcount = 1,outxdim ! Loop over x and y
    Do ycount = 1,outydim
      If (outdata(xcount,ycount).lt.laps_missing) Then
        num_good = num_good + 1
        If (outdata(xcount,ycount).gt.data_max) data_max = outdata(xcount,ycount)
        If (outdata(xcount,ycount).lt.data_min) data_min = outdata(xcount,ycount)
        data_mean = ((data_mean * (num_good-1)) + outdata(xcount,ycount)) / num_good
      EndIf
    EndDo
    EndDo ! Loop over x and y
    
    If (debug.gt.1) Then ! Output grid information
      write(6,'(a)') ''
      write(6,'(a)') 'Output information '
      write(6,'(a)') '^^^^^^^^^^^^^^^^^^ '
      write(6,'(a,i4,a,i4)') 'Dimensions : ',outxdim,' x ',outydim
      write(6,'(a,f7.2,a,f7.2,a,f7.2,a,f7.2)') 'Lat Corners: ',outlat(1,1),', ', &
           outlat(outxdim,1),', ',outlat(1,outydim),', ',outlat(outxdim,outydim)
      write(6,'(a,f7.2,a,f7.2,a,f7.2,a,f7.2)') 'Lat Corners: ',outlon(1,1),', ', &
           outlon(outxdim,1),', ',outlon(1,outydim),', ',outlon(outxdim,outydim)
      write(6,'(a,i7,a,i7,a,f6.2,a)') 'Number of good values/total (%) : ', num_good, &
              ' / ',outxdim*outydim,'  (', (num_good*100.0)/(outxdim*outydim),'% )'
      write(6,'(a,f6.3,a,f6.3,a,f7.3)') 'Data in mm (w/o missing) min/mean/max : ', &
                       1000.0*data_min,' / ',1000.0*data_mean,' / ',1000.0*data_max
    EndIf

    If (data_mean.eq.0.0) write(6,'(a)') ''
    If (data_mean.eq.0.0) write(6,'(a)') '***WARNING***'
    If (data_mean.eq.0.0) write(6,'(a)') '   The output data mean is zero! Check that your'
    If (data_mean.eq.0.0) write(6,'(a)') '   Check that your LAPS_DATA_ROOT is set up correctly'
    If (data_mean.eq.0.0) write(6,'(a)') '***WARNING***'
    If (data_mean.eq.0.0) write(6,'(a)') ''

    ! Output the files using the LAPS standards
    directory = './'
    extension = 'st4'
    variable  = 'PPT'
    units     = 'mm'
    comment   = 'Stage IV precipitation data (mm) remapped to LAPS domain'
    CALL write_laps_data(lapsi4t,directory,extension,outxdim,outydim,1,1, &
                         variable,lvl,lvl_coord,units,comment,outdata,err)
    err = 0
    RETURN
  END SUBROUTINE st4_driver

  ! The following are used by the st4_driver

  !==============================================================================================!
  ! Subroutine cdf(status, err)                                                                  !
  !==============================================================================================!
  !---   Purpose - Error handling framework for netCDF                                           !
  !----------------------------------------------------------------------------------------------!
  SUBROUTINE cdf(status, err)
    USE netcdf
    Integer, intent (In)  :: status
    Integer, intent (Out) :: err
      err = 0
    If (status.ne.nf90_noerr) Then 
      err = 1
      write(6,'(a)') 'From netCDF: '//nf90_strerror(status)
      RETURN
    EndIf
  END SUBROUTINE cdf  

  !==============================================================================================!
  !  SUBROUTINE | alloc_error(errcode)                                                           !
  !==============================================================================================!
  !---   Purpose - Displays information about a memory allocation error                          !
  !---   Input - error code                                                                      !
  !----------------------------------------------------------------------------------------------!
  SUBROUTINE alloc_error(errcode)
    Integer, Intent(In) :: errcode
    Write(6,'(a)') 'An error occurred while allocating memory.'
    if (errcode.eq.1) Write(6,'(a)') 'Error in system routine attempting to do allocation'
    if (errcode.eq.2) Write(6,'(a)') 'An invalid data object has been specified for allocation'
    if (errcode.eq.3) Write(6,'(a)') 'Error in system routine attempting to do allocation'
    if (errcode.eq.4) Write(6,'(a)') 'An invalid data object has been specified for allocation'
    STOP
  END SUBROUTINE alloc_error
