!  output in netCDF format (assumes 1 forecast time per
!  file)
!
!  Brent L. Shaw, Weathernews Americas, Inc., 2005

SUBROUTINE lga_driver_wrfarw(nx_laps,ny_laps,nz_laps,bgpath, cmodel, & ! add nx_laps,ny_laps,nz_laps by Wei-Ting(130326)
                             use_analysis,forecast_length, &
                             luse_sfc_bkgd, &
                             i4time_now, smooth_fields,lga_status)

  USE wrf_netcdf
  USE time_utils
  USE wrf_lga
  IMPLICIT NONE
  include 'netcdf.inc'
  ! === add variables by Wei-Ting(130326) ===
  INTEGER, INTENT(IN)          :: nx_laps,ny_laps,nz_laps ! LAPS grid dimensions
  REAL                         :: ht_1d(nz_laps), &       !LAPS vert grid (SIGMA_HT grid only)
                                  pr1d_pa(nz_laps), &     !LAPS pressures (pa)
                                  pr1d_mb(nz_laps)        !LAPS pressures (mb)
  CHARACTER(LEN=9)             :: wfo_fname13_to_fname9
  CHARACTER(LEN=9)             :: fname9_reftime          ! reftime julian date
  INTEGER                      :: i4time_reftime,i4time_fcst(2)
  CHARACTER(LEN=256)           :: outdir
  CHARACTER(LEN=31)            :: ext
  INTEGER                      :: iloop,len_dir,hr_now,min_now
  ! === end of add by Wei-Ting(130326) ===
  CHARACTER(LEN=*),INTENT(IN)  :: bgpath
  CHARACTER(LEN=12)            :: cmodel
  LOGICAL, INTENT(IN)          :: use_analysis
  INTEGER, INTENT(IN)          :: forecast_length ! hours
  LOGICAL, INTENT(IN)          :: luse_sfc_bkgd
  INTEGER, INTENT(IN)          :: i4time_now 
  LOGICAL, INTENT(IN)          :: smooth_fields
  INTEGER, INTENT(OUT)         :: lga_status ! 1 = success

  ! Local variables
  INTEGER              :: nfiles,cdf,istatus
  INTEGER              :: files_i4time(2)  
  CHARACTER(LEN=256)   :: filenames(3) ! plus 1 dim. by Wei-Ting (130312) to put previous time
  CHARACTER(LEN=19)    :: reftime
  REAL                 :: dt
  INTEGER              :: itimestep, tau_hr, tau_min, tau_sec 
  LOGICAL              :: already_done
  lga_status = 1  ! If we encounter a failure, change to 0


  ! Get list of acceptable files
  CALL get_acceptable_wrf(bgpath,i4time_now,nfiles, &
     files_i4time, filenames)


  ! Do we have an exact time match?  
  IF (nfiles .EQ. 1) THEN
    ! Yes, get its attributes
    PRINT *, "Found exact match: ", TRIM(filenames(1))
    ! Make sure file is ready
    CALL wrfio_wait(filenames(1),300)
    CALL open_wrfnc(filenames(1),cdf,istatus)
    IF (istatus .NE. 0) THEN
      PRINT *, "Problem opening ",TRIM(filenames(1))
      lga_status = 0
      RETURN 
    ENDIF
    ! Get time information
    CALL get_wrf2_timeinfo(cdf,reftime,dt,itimestep,tau_hr, &
           tau_min, tau_sec, istatus)
    CALL close_wrfnc(cdf)
    IF (istatus .NE. 0) THEN
      PRINT*, "Problem getting time info from ", TRIM(filenames(1))
      lga_status = 0
      RETURN
    ENDIF
    ! Is it a new enough forecast?
    IF (tau_hr .LE. forecast_length) THEN
      PRINT *, "Forecast hour = ",tau_hr

      ! Only process if we don't already have an LGA output file from WRFARW
      ! that matches this reftime+tau
      PRINT *, "Calling check_wrf_lga for:",cmodel,reftime
      CALL check_wrf_lga(cmodel,reftime,tau_hr,tau_min,already_done) 
      print *, "Already done: ",already_done
      IF (.NOT. already_done) THEN
        PRINT *, "Calling wrf2lga"
        CALL wrf2lga(filenames(1:2),i4time_now,cmodel,istatus) ! modified by Wei-Ting (130326) to insert previous time ( filename(1) -> filname(1:2) )
        IF (istatus .NE. 1) THEN
          PRINT *, "Failure in wrf2lga"
          lga_status = 0
          RETURN
        ENDIF
      ELSE
        PRINT *, "File already processed previously"
        lga_status = 1
        RETURN
      ENDIF
    ELSE
      PRINT *, "Forecast is too old :", tau_hr, forecast_length
      lga_status = 0
      RETURN
    ENDIF
  ELSEIF(nfiles .EQ. 2) THEN
    ! No, but we do have 2 files that bound this time
    PRINT *, "Found bounding files: "
    PRINT *, TRIM(filenames(2)) ! eariler time : exchange two lines by Wei-Ting(130326)
    PRINT *, TRIM(filenames(1)) ! later time   : exchange two lines by Wei-Ting(130326)
    PRINT *, " "
    ! Does latest bound exceed forecast limit?
             ! No...read both in, then time interpolate
                ! Get domain info
                ! Allocate three sets of arrays (t-1,t,t+1)
                ! Get destaggered variables for t-1 and t-1
                ! Time interpolate to t
                ! Deallocate t-1 and t+1

             ! Yes...exit with error
        ! No
          ! Exit with error
    
    ! ===== Below codes is added to complete        =====
    ! ===== time interpolation by Wei-Ting (130326) =====
    DO iloop = 2,1,-1
       PRINT *, "Processing :", TRIM(filenames(iloop))
       ! Make sure file is ready
       CALL wrfio_wait(filenames(iloop),300)
       CALL open_wrfnc(filenames(iloop),cdf,istatus)
       IF (istatus .NE. 0) THEN
         PRINT *, "Problem opening ",TRIM(filenames(iloop))
         lga_status = 0
         RETURN 
       ENDIF
       ! Get time information
       CALL get_wrf2_timeinfo(cdf,reftime,dt,itimestep,tau_hr, &
              tau_min, tau_sec, istatus)
       CALL close_wrfnc(cdf)
       IF (istatus .NE. 0) THEN
         PRINT*, "Problem getting time info from ", TRIM(filenames(iloop))
         lga_status = 0
         RETURN
       ENDIF
       ! Is it a new enough forecast?
       IF (tau_hr .LE. forecast_length) THEN
         PRINT *, "Forecast hour = ",tau_hr

         ! Only process if we don't already have an LGA output file from WRFARW
         ! that matches this reftime+tau
         PRINT *, "Calling check_wrf_lga for:",cmodel,reftime
         CALL check_wrf_lga(cmodel,reftime,tau_hr,tau_min,already_done) 
         print *, "Already done: ",already_done
         IF (.NOT. already_done) THEN
           PRINT *, "Calling wrf2lga"
           CALL wrf2lga(filenames(iloop:iloop+1),i4time_now,cmodel,istatus)
           IF (istatus .NE. 1) THEN
             PRINT *, "Failure in wrf2lga"
             lga_status = 0
             RETURN
           ENDIF
         ELSE
           PRINT *, "File already processed previously"
           lga_status = 1
           ! RETURN
           CONTINUE
         ENDIF
       ELSE
         PRINT *, "Forecast is too old :", tau_hr, forecast_length
         lga_status = 0
         RETURN
       ENDIF
    ENDDO
    ! ===== Time Interpolation =====
    fname9_reftime = wfo_fname13_to_fname9&
                  (reftime(1:4)//reftime(6:7)//reftime(9:13)//reftime(15:16))
    call i4time_fname_lp(fname9_reftime,i4time_reftime,istatus)
    i4time_fcst(1) = files_i4time(1)-i4time_reftime
    i4time_fcst(2) = files_i4time(2)-i4time_reftime
    hr_now = (i4time_now-i4time_reftime)/3600
    min_now = MOD(i4time_now-i4time_reftime,3600)/60
    PRINT *, "Calling check_wrf_lga for:",cmodel,reftime
    CALL check_wrf_lga(cmodel,reftime,hr_now,min_now,already_done)
    print *, "Already done: ",already_done
    IF (.NOT. already_done) THEN
       print*,'get 1d pressures'
       call get_pres_1d(i4time_now,nz_laps,pr1d_pa,istatus)
       if(istatus.ne.1)then
          print*,'Error returned from get_pres_1d'
          print*,'Check pressures.nl or nk_laps in nest7grid.parms'
          stop
       endif
       pr1d_mb(:)=pr1d_pa(:)/100.  ! Pa to mb
       print*,i4time_now,i4time_reftime,i4time_reftime, &
              i4time_fcst(2),i4time_fcst(1), &
              files_i4time(2),files_i4time(1)
       ext = 'lga'
       call get_directory(ext,outdir,len_dir)
       print*,outdir,ext,nz_laps
       ! interp 3D fields
       call time_interp(outdir,ext, &
                   nx_laps,ny_laps,nz_laps,6,pr1d_mb,ht_1d, &
                   files_i4time(2),files_i4time(1), &
                   i4time_now,i4time_reftime,i4time_fcst(1), &
                   i4time_reftime,i4time_fcst(2))
       ext = 'lgb'
       call get_directory(ext,outdir,len_dir)
       print*,outdir,ext
       ! interp 2D fields
       call time_interp(outdir,ext, &
                   nx_laps,ny_laps,1,10,pr1d_mb(1),ht_1d(1), &
                   files_i4time(2),files_i4time(1), &
                   i4time_now,i4time_reftime,i4time_fcst(1), &
                   i4time_reftime,i4time_fcst(2))
       
    ELSE
       PRINT *, "File already processed previously"
       lga_status = 1
       RETURN
    ENDIF
    ! ===== End of modification by Wei-Ting (130326) =====
  ELSE
    PRINT *, "No acceptable WRF files found!"
    lga_status = 0
    RETURN
  ENDIF

  RETURN
END SUBROUTINE lga_driver_wrfarw

SUBROUTINE  check_wrf_lga(cmodel,reftime,tau_hr,tau_min,already_done)
  IMPLICIT NONE
  include 'netcdf.inc'
  CHARACTER(LEN=12),INTENT(IN)         :: cmodel
  CHARACTER(LEN=19), INTENT(IN)        :: reftime
  INTEGER,INTENT(IN)                   :: tau_hr,tau_min
  LOGICAL,INTENT(OUT)                  :: already_done

  CHARACTER(LEN=9),EXTERNAL            :: wfo_fname13_to_fname9
  CHARACTER(LEN=9)                     :: lga_reftime
  CHARACTER(LEN=13)                    :: reftime13
  CHARACTER(LEN=200)                   :: lgadir
  CHARACTER(LEN=13)                    :: lgatime 
  CHARACTER(LEN=200)                   :: lgafile
  CHARACTER(LEN=132)                   :: ht_comment
  INTEGER                              :: lendir,cdf,nstatus,varid,cmodel_len
  INTEGER                              :: start(3), count(3)
  already_done = .false.
  PRINT *, "In check_wrf_lga:"
  PRINT *, "cmodel = ", cmodel
  PRINT *, "reftime = ", reftime
  PRINT *, "Tau_hr / Tau_min = " , tau_hr, tau_min   
  reftime13 = reftime(1:4)//reftime(6:7)//reftime(9:13)//reftime(15:16)
  print *, reftime13
  lga_reftime = wfo_fname13_to_fname9(reftime13)
  print *, reftime13,lga_reftime
  WRITE(lgatime, '(A9,I2.2,I2.2)') lga_reftime,tau_hr,tau_min
  CALL get_directory('lga',lgadir,lendir)
  lgafile = lgadir(1:lendir)//lgatime//'.lga'
  print *, "lgafile ", lgafile
  INQUIRE(FILE=TRIM(lgafile),EXIST=already_done)
  IF (already_done) THEN

    ! Check the comment field to see if this is WRFARW
    nstatus = NF_OPEN(lgafile,NF_NOWRITE,cdf)
    IF (nstatus .NE. NF_NOERR) THEN
      already_done = .FALSE.
      RETURN
    ENDIF
    nstatus = NF_INQ_VARID(cdf, 'ht_comment', varid)
    IF (nstatus .NE. NF_NOERR) THEN
      PRINT *, "Problem getting ht_comment"
      nstatus = NF_CLOSE(cdf)
      already_done = .false.
      RETURN
    ENDIF
    start = (/1,1,1/)
    count = (/132,1,1/)
    nstatus = NF_GET_VARA_TEXT(cdf,varid,start,count,ht_comment)
    nstatus = NF_CLOSE(cdf)
    IF (cmodel(1:6) .EQ. ht_comment(1:6)) THEN
      PRINT *, "File already processed for ",cmodel
      already_done = .true.
    ELSE
      print *, cmodel(1:6),"|",ht_comment(1:6)
      PRINT *, "Different cmodel : ", ht_comment(1:6) 
      already_done = .FALSE.
    ENDIF
  ENDIF
  RETURN
END SUBROUTINE check_wrf_lga
