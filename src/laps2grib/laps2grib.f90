!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Program: 	laps2grib
!!
!! Purpose: 	Converts a subset of LAPS analysis grids to a grib2 file
!!
!! Author:	Brent Shaw, Weathernews Inc.
!!
!! History:
!!		07 Dec 2006	Brent Shaw	Initial version
!!		25 Nov 2011 	Paula McCaslin	Modified to add accumulated varaibles
!!		15 Dec 2011 	Paula McCaslin	Modified to add model varaibles (from fsf, fua)
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM laps2grib

  !! Requires some modules in LAPS_SRC_ROOT/lib/modules
  USE grib2
  !! Local modules
  USE laps2grib_config
  USE lapsdata
  include 'lapsparms.cmn'     ! laps_cycle_time, model_cycle_time, model_fcst_intvl
  
  IMPLICIT NONE

  INTEGER                     :: istatus,year,month,day,hour,minute,second
  INTEGER                     :: filestatus, g2lun
  INTEGER                     :: ncid,varid,ncstat
  INTEGER                     :: k,n
  INTEGER, PARAMETER          :: data_type = 0
  INTEGER                     :: reftime_sig !0=analysis, 1=start of fcst
  INTEGER, PARAMETER          :: process_type = 0
  INTEGER, PARAMETER          :: bg_process_id = 255
  INTEGER, PARAMETER          :: cutoff_hr = 0
  INTEGER, PARAMETER          :: cutoff_min =0 
  INTEGER                     :: time_unit_indicator
  INTEGER                     :: time_range
  INTEGER, PARAMETER          :: levtype_iso = 100 !fixed level type
  INTEGER, PARAMETER          :: levscale_iso = 0
  INTEGER, PARAMETER          :: g2miss = 255
  INTEGER, PARAMETER          :: newrec = 1 
  INTEGER, PARAMETER          :: pack_method = 2
  INTEGER                     :: miss_mgmt,inomiss
  INTEGER                     :: val_time, ref_time, accum_time
  INTEGER                     :: eyear,emonth,eday,ehour,eminute,esecond
  INTEGER                     :: etime_value, etime_unit
  
  REAL, ALLOCATABLE           :: lapsdata2d(:,:),lapsdata3d(:,:,:),slab(:,:)
  REAL                        :: r_missing, rhhmm
  CHARACTER(LEN=32),PARAMETER :: vtag = "laps2grib V1.2, 01 Dec 2011, WNI" !"laps2grib V1.0, 08 Dec 2006, WNI" 
  CHARACTER(LEN=256)          :: laps_data_root,lapsfile
  CHARACTER(LEN=256)          :: g2file, g2file_tmp
  CHARACTER(LEN=512)          :: syscmd
 
  INTEGER                     :: i, iargc, max_args, ihhmm
  INTEGER                     :: modeltime_passed_in
  CHARACTER(LEN=100)          :: vtab, forecast_id
  CHARACTER(LEN=5)            :: hhmm
  CHARACTER(LEN=14)           :: file_a9time, modeltime
  LOGICAL                     :: dir_exists

  ! Print banner
  print *, "=================================================================="
  print *, "********** ",vtag," **********"
  print *, ""
  print *, " USAGE:       laps2grib.exe [vtab]"
  print *, " MODEL USAGE: laps2grib.exe vtab [hh]hmm forecast_id [modeltime]"
  print *, "        (e.g. laps2grib.exe wfr2grib.vtab 1200 wrf-hrrr 130451800)"
  print *, "        (optional modeltime has the format YYJJJHHMM)"
  print *, "=================================================================="

  vtab = 'laps2grib.vtab'
  reftime_sig = 0
  forecast_id = ''
  hhmm = ''
  max_args = iargc()

  IF (max_args .EQ. 1) THEN
    ! expect vtable name
    CALL GETARG(1,vtab)

  ELSE IF ((max_args .EQ. 3) .OR. (max_args .EQ. 4)) THEN
    reftime_sig = 1
    ! expect vtable name
    CALL GETARG(1,vtab)

    ! expect hhmm for forecast time
    CALL GETARG(2,hhmm)
    READ(hhmm,*,IOSTAT=ihhmm) rhhmm 
    IF ( rhhmm .EQ. 0. .AND. ihhmm .NE. 0) THEN
        PRINT *, "The usage requires HHMM for the forecast time, this is not an expected number: ", hhmm
        STOP 
    ELSE IF (rhhmm .EQ. 0.) THEN 
        hhmm = '0000'
    ELSE IF ( LEN_TRIM(hhmm) .LT. 3) THEN
        PRINT *, "The usage requires HHMM for the forecast time, the string is too short: ", hhmm
        STOP 
    ELSE IF ( LEN_TRIM(hhmm) .LE. 3 .AND. index(hhmm,' ') .NE. 0) THEN
        PRINT *, "--> Index Value", index(hhmm,' ')
        hhmm = '0'//hhmm
    ENDIF

    ! expect ensemble forecast_id name, .e.g. mean, or wrf-hrrr
    CALL GETARG(3,forecast_id)
    forecast_id = '/'//forecast_id
    modeltime_passed_in = 0
    IF (max_args .EQ. 4) THEN
      CALL GETARG(4,modeltime)
      IF ( LEN_TRIM(modeltime) .NE. 9) THEN
        print *, 'Modeltime passed in is not the correct length: ',trim(modeltime)
        print *, 'Defaulting to time in file modeltime.dat'
      ELSE
        modeltime_passed_in = 1
      ENDIF
    ENDIF
  ELSE IF (max_args .NE. 0) THEN
        STOP "Check the usage statement above..."
  ENDIF

  ! Get the LAPS_DATA_ROOT
  CALL GETENV('LAPS_DATA_ROOT',laps_data_root)
  IF ( LEN_TRIM(laps_data_root) .LT. 1) THEN
    PRINT *, "LAPS_DATA_ROOT not set!"
    STOP 'no LAPS_DATA_ROOT'
  ENDIF

  PRINT *, "-- LAPS_DATA_ROOT=",TRIM(laps_data_root)
 
  ! Read the laps2grib namelist
  CALL read_laps2grib_nl(laps_data_root)

  if(.not. lrun_laps2grib)then
      write(6,*)'lrun_laps2grib is set to FALSE - stopping program'
      stop
  endif

  ! Set up projection
  CALL get_laps_proj(laps_data_root) 
  CALL get_laps_plevels(laps_data_root)

  ! Get laps analysis valid time
  IF (LEN_TRIM(hhmm) .EQ. 0) THEN 
     CALL get_laps_analtime
  ELSE
     CALL get_laps_modeltime(modeltime,modeltime_passed_in)
  ENDIF
  CALL cv_i4tim_int_lp(i4time,year,month,day,hour,minute,second,istatus)
  year = year + 1900
  file_a9time=a9time
  IF (hhmm .NE. '') file_a9time=a9time//hhmm
  PRINT *, "-- Using Timestamp: ", file_a9time

  ! Configure the variable list
  CALL get_data_config(laps_data_root,vtab)

  ! Get the rmissing value
  CALL get_r_missing_data(r_missing,istatus)
  PRINT *, "-- R_MISSING: ",r_missing

  ! Check the GRIB dir
  g2file = TRIM(output_path)//TRIM(forecast_id)
  INQUIRE(FILE=g2file,EXIST=dir_exists)
  IF (.NOT. dir_exists) THEN
      PRINT *, "Dir does not exist: ",TRIM(g2file)
      STOP 
  ENDIF

  ! Open the GRIB file
  g2file = TRIM(output_path)//TRIM(forecast_id)//'/'//TRIM(file_a9time)//'.gr2'
  g2file_tmp = TRIM(output_path)//TRIM(forecast_id)//'/.'//TRIM(file_a9time)//'.gr2'
  CALL init_grib2_file(g2file_tmp,laps_proj,center_id,subcenter_id, &
                       reftime_sig,year,month,day,hour,minute,second, &
                       prod_status,data_type,g2lun,istatus)

  ! Process any 3D isobaric variables
  IF (n_iso .GT. 0) THEN
    ALLOCATE(lapsdata3d(nx,ny,nz))
    ALLOCATE(slab(nx,ny))
    PRINT *, "-- Processing 3D Isobaric Fields"
    loop3d: DO n = 1, n_iso

      IF (meta_iso3d(n)%qbal .EQ. 0) THEN
        lapsfile = TRIM(laps_data_root)//'/lapsprd/'//meta_iso3d(n)%ext//TRIM(forecast_id)//'/'//TRIM(file_a9time)//'.'// &
                         meta_iso3d(n)%ext
      ELSE
        lapsfile = TRIM(laps_data_root)//'/lapsprd/balance/'//meta_iso3d(n)%ext//'/'//&
                     TRIM(file_a9time)//'.'//meta_iso3d(n)%ext
      ENDIF

      filestatus = 0
      CALL ncread_3d(lapsfile,nx,ny,nz,meta_iso3d(n)%var,lapsdata3d,val_time,ref_time,istatus)
      IF (istatus .EQ. 1) THEN
        filestatus = 1  ! ncread_3d returned data...OK to write out grib2
        CALL calc_val_ref_times(val_time,ref_time,time_unit_indicator,time_range,istatus)
        IF (istatus .NE. 1) THEN
          filestatus = 0  ! could not calc val_ref times...DO NOT write out grib2
          PRINT *, "Could not calc val_ref times: ",meta_iso3d(n)%ext," VAR ",meta_iso3d(n)%var
        ENDIF
      ELSE
        PRINT *, "Could not read file: ",meta_iso3d(n)%ext," VAR ",meta_iso3d(n)%var
      ENDIF

      IF (filestatus .EQ. 1) THEN
        press_loop:  DO k = 1, nz
          IF( (plevels_pa(k) .GE. meta_iso3d(n)%pmin).AND. &
              (plevels_pa(k) .LE. meta_iso3d(n)%pmax) ) THEN
         
            slab = lapsdata3d(:,:,k)
            IF ( MAXVAL(slab) .EQ. r_missing ) THEN
              inomiss = 0
              miss_mgmt = 1
              IF (meta_iso3d(n)%conv_fac .NE. 1.) THEN
                WHERE(slab .NE. r_missing) slab = slab * meta_iso3d(n)%conv_fac
              ENDIF

            ELSE
              inomiss = 1
              miss_mgmt = 0
              IF (meta_iso3d(n)%conv_fac .NE. 1.) THEN
                slab = slab * meta_iso3d(n)%conv_fac
              ENDIF
            ENDIF
            PRINT '("GRIB: ",A3,1x,A3,1x,I5,"mb",I3,I3,I4,F12.4,F12.4)',meta_iso3d(n)%ext, &
                meta_iso3d(n)%var,NINT(plevels_pa(k)*0.01),&
                meta_iso3d(n)%discipline,meta_iso3d(n)%cat_id,meta_iso3d(n)%parm_id,&
                MINVAL(slab),MAXVAL(slab)
            CALL write_grib2_template0(g2lun,meta_iso3d(n)%discipline, &
                 meta_iso3d(n)%cat_id,meta_iso3d(n)%parm_id, process_type, &
                 bg_process_id,process_id,cutoff_hr,cutoff_min,time_unit_indicator, &
                 time_range, levtype_iso,levscale_iso, NINT(plevels_pa(k)),g2miss,g2miss,&
                 g2miss,pack_method,meta_iso3d(n)%scale_fac,miss_mgmt, &
                 nx,ny,newrec,inomiss,r_missing, r_missing,slab)
          ELSE  
            PRINT '("SKIP: ",A3,1x,A3,1x,I5,"mb")', meta_iso3d(n)%ext,meta_iso3d(n)%var, &
                   NINT(plevels_pa(k)*0.01)
          ENDIF
        ENDDO press_loop

        
      ELSE  
        PRINT *, "Problem getting: ",TRIM(forecast_id)," EXT .",meta_iso3d(n)%ext," VAR ",meta_iso3d(n)%var
      ENDIF
    ENDDO loop3d
    DEALLOCATE(lapsdata3d)
    DEALLOCATE(slab)
  ENDIF

  ! Process any 2D variables
  IF (n_2d .GT. 0) THEN
    ALLOCATE(lapsdata2d(nx,ny))
    PRINT *, "-- Processing 2D Fields"
    loop2d: DO n =1, n_2d
      IF (meta_2d(n)%ext .NE. 'n7g') THEN
        lapsfile = TRIM(laps_data_root)//'/lapsprd/'//meta_2d(n)%ext//TRIM(forecast_id)//'/'//TRIM(file_a9time)//'.'// &
                         meta_2d(n)%ext

      ELSE
        lapsfile = TRIM(laps_data_root)//'/static/static.nest7grid'
      ENDIF

      filestatus = 0
      CALL ncread_2d(lapsfile,nx,ny,meta_2d(n)%var,lapsdata2d,val_time,ref_time,istatus)
      IF (istatus .EQ. 1) THEN
        filestatus = 1  ! ncread_2d returned data...OK to write out grib2
        CALL calc_val_ref_times(val_time,ref_time,time_unit_indicator,time_range,istatus)
        IF (istatus .NE. 1) THEN
          filestatus = 0  ! could not calc val_ref times...DO NOT write out grib2
          PRINT *, "Could not calc val_ref times: ",meta_2d(n)%ext," VAR ",meta_2d(n)%var
        ENDIF
      ELSE
        PRINT *, "Could not read file: ",meta_2d(n)%ext," VAR ",meta_2d(n)%var
      ENDIF

      IF (filestatus .EQ. 1) THEN
        IF ( MAXVAL(lapsdata2d) .EQ. r_missing ) THEN
          inomiss = 0
          miss_mgmt = 1
          IF (meta_2d(n)%conv_fac .NE. 1.) THEN
            WHERE(lapsdata2d .NE. r_missing) lapsdata2d = lapsdata2d * meta_2d(n)%conv_fac
          ENDIF
                                                                                                  
        ELSE
          inomiss = 1
          miss_mgmt = 0
          IF (meta_2d(n)%conv_fac .NE. 1.) THEN
            lapsdata2d = lapsdata2d * meta_2d(n)%conv_fac
          ENDIF
        ENDIF

        PRINT '("GRIB: ",A3,1x,A3,1x,F12.4,F12.4)',&
          meta_2d(n)%ext, meta_2d(n)%var, MINVAL(lapsdata2d),MAXVAL(lapsdata2d)
        CALL write_grib2_template0(g2lun,meta_2d(n)%discipline, &
                 meta_2d(n)%cat_id,meta_2d(n)%parm_id, process_type, &
                 bg_process_id,process_id,cutoff_hr,cutoff_min,time_unit_indicator, &
                 time_range, meta_2d(n)%lev1_type,meta_2d(n)%lev1_scale,meta_2d(n)%lev1_value,&
                 meta_2d(n)%lev2_type,meta_2d(n)%lev2_scale,meta_2d(n)%lev2_value,&
                 pack_method,meta_2d(n)%scale_fac,miss_mgmt, &
                 nx,ny,newrec,inomiss,r_missing, r_missing,lapsdata2d)

      ELSE 
        PRINT *, "Problem getting: ",TRIM(forecast_id)," EXT .",meta_2d(n)%ext," VAR ",meta_2d(n)%var
      ENDIF
    ENDDO loop2d
    DEALLOCATE(lapsdata2d)
  ENDIF

  ! Process any Accumulated 2D variables
  IF (n_accum2d .GT. 0) THEN
    ALLOCATE(lapsdata2d(nx,ny))
    PRINT *, "-- Processing Accumulated 2D Fields"
    loop_a2d: DO n =1, n_accum2d
      lapsfile = TRIM(laps_data_root)//'/lapsprd/'//meta_accum2d(n)%ext//TRIM(forecast_id)//'/'//TRIM(file_a9time)//'.'// &
                         meta_accum2d(n)%ext

      filestatus = 0
      CALL ncread_2d(lapsfile,nx,ny,meta_accum2d(n)%var,lapsdata2d,val_time,ref_time,istatus)
      PRINT *, "ncread_2d(", lapsfile,nx,ny,n, val_time,ref_time,istatus, ")"
      IF (istatus .EQ. 1) THEN
        filestatus = 1  ! ncread_2d returned data...OK to write out grib2
        CALL calc_val_ref_times(val_time,ref_time,time_unit_indicator,time_range,istatus)
        IF (istatus .NE. 1) THEN
          filestatus = 0  ! could not calc val_ref times...DO NOT write out grib2
          PRINT *, "Could not calc val_ref times: ",meta_accum2d(n)%ext," VAR ",meta_accum2d(n)%var
        ENDIF
      ELSE
        PRINT *, "Could not read file: ",meta_accum2d(n)%ext," VAR ",meta_accum2d(n)%var
      ENDIF

      IF (filestatus .EQ. 1) THEN
        IF ( MAXVAL(lapsdata2d) .EQ. r_missing ) THEN
          inomiss = 0
          miss_mgmt = 1
          IF (meta_accum2d(n)%conv_fac .NE. 1.) THEN
            WHERE(lapsdata2d .NE. r_missing) lapsdata2d = lapsdata2d * meta_accum2d(n)%conv_fac
          ENDIF
                                                                                                  
        ELSE
          inomiss = 1
          miss_mgmt = 0
          IF (meta_accum2d(n)%conv_fac .NE. 1.) THEN
            lapsdata2d = lapsdata2d * meta_accum2d(n)%conv_fac
          ENDIF
        ENDIF
        
        CALL calc_accum_time(lapsfile,meta_accum2d(n)%ext,meta_accum2d(n)%var,val_time,ref_time,accum_time,&
			   time_unit_indicator,time_range,etime_unit,etime_value,istatus)
        IF (istatus .EQ. 0) PRINT *, "Problem with CALC_ACCUM_TIME"

        PRINT *, 'after calc_accum_time ', ref_time, ', ', ref_time, ' ', &
           accum_time
        PRINT *, 'etime: ', etime_unit, ', ', etime_value

        ! Octet 49-50, set to current runtime
        meta_accum2d(n)%etime_unit = etime_unit 
        meta_accum2d(n)%etime_value= etime_value

        ! cv_i4tim_int_lp expects LAPS time  KLH - 20 MAR 2014
        val_time = val_time + 315619200 
	! --- Valid date, time to describe data for GRIB2.
        CALL cv_i4tim_int_lp(val_time,eyear,emonth,eday,ehour,eminute,esecond,istatus)
        !PRINT *, 'cv_i4tim_int_lp(', val_time,eyear,emonth,eday,ehour,&
        !eminute,esecond,istatus, ')'
        ! Add 10 to eyear to have time in grib2 files to account for 10 yr difference between i4time and unixtime
        ! Octet 35-40, set to current runtime
        !meta_accum2d(n)%eyear = eyear + 1900 + 10
        ! This 10 year offset is not needed when valid LAPS time is passed to cv_itim_int_lp   KLH -- 20 MAR 2014

        meta_accum2d(n)%eyear = eyear + 1900
        meta_accum2d(n)%emon  = emonth
        meta_accum2d(n)%eday  = eday
        meta_accum2d(n)%ehour = ehour
        meta_accum2d(n)%emin  = 0

        PRINT '("GRIB: ",A3,1x,A3,1x,F12.4,F12.4)',&
          meta_accum2d(n)%ext, meta_accum2d(n)%var, MINVAL(lapsdata2d),MAXVAL(lapsdata2d)
        CALL write_grib2_template8(g2lun,meta_accum2d(n)%discipline, &
                 meta_accum2d(n)%cat_id,meta_accum2d(n)%parm_id, process_type, &
                 bg_process_id,process_id,cutoff_hr,cutoff_min,time_unit_indicator, &
                 time_range, meta_accum2d(n)%lev1_type,meta_accum2d(n)%lev1_scale,meta_accum2d(n)%lev1_value,&
                 meta_accum2d(n)%lev2_type,meta_accum2d(n)%lev2_scale,meta_accum2d(n)%lev2_value,&
                 meta_accum2d(n)%eyear,meta_accum2d(n)%emon,meta_accum2d(n)%eday,&
                 meta_accum2d(n)%ehour,meta_accum2d(n)%emin,meta_accum2d(n)%esec,&
                 meta_accum2d(n)%ntimes,meta_accum2d(n)%ntimes_miss,&
                 meta_accum2d(n)%stattype,meta_accum2d(n)%periodtype,&
                 meta_accum2d(n)%etime_unit,meta_accum2d(n)%etime_value,&
                 pack_method,meta_accum2d(n)%scale_fac,miss_mgmt,&
                 nx,ny,newrec,inomiss,r_missing, r_missing,lapsdata2d)
      ELSE 
        PRINT *, "Problem getting: ",TRIM(forecast_id)," EXT .",meta_accum2d(n)%ext," VAR ",meta_accum2d(n)%var
      ENDIF
    ENDDO loop_a2d
    DEALLOCATE(lapsdata2d)
  ENDIF

  ! Close the file and rename it to the final output name
  PRINT *, "-- Closing File"
  CALL close_grib2_file(g2lun)
  syscmd = 'mv '// TRIM(g2file_tmp) // ' ' // TRIM(g2file)   
  PRINT *, TRIM(syscmd)
  CALL system(syscmd)
  PRINT *, "======================================================"
  PRINT *, "**********      laps2grib completed         **********"
  PRINT *, "======================================================"
END PROGRAM laps2grib

!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ncread_3d(ncfile,nx,ny,nz,ncvar,data3d,nc_val,nc_ref,istatus)

    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN)  :: ncfile
    INTEGER, INTENT(IN)           :: nx,ny,nz
    CHARACTER(LEN=3), INTENT(IN)  :: ncvar
    REAL,INTENT(OUT)              :: data3d(nx,ny,nz)
    INTEGER, INTENT(OUT)          :: istatus

    INTEGER  :: ncid,vid,ncstat
    INTEGER  :: nc_val,nc_ref
    LOGICAL  :: file_exists
    INCLUDE 'netcdf.inc'

    istatus = 1
    INQUIRE(FILE=ncfile,EXIST=file_exists)
    IF (.NOT. file_exists) THEN
      PRINT *, "File not found: ",TRIM(ncfile)
      istatus =0
      RETURN
    ENDIF
    ncstat = NF_OPEN(ncfile,NF_NOWRITE,ncid)
    IF (ncstat .NE. NF_NOERR) THEN
      PRINT *, "Error opening: ",TRIM(ncfile)
      istatus = 0
      RETURN
    ENDIF
    ncstat = NF_INQ_VARID(ncid,ncvar,vid)
    IF (ncstat .EQ. NF_NOERR) THEN
      ncstat = NF_GET_VAR_REAL(ncid,vid,data3d)
      IF (ncstat .NE. NF_NOERR) THEN
        PRINT *, "Problem getting data for ", TRIM(ncvar)
        istatus = 0
      ENDIF 
    ELSE
      PRINT *,"Could not find Var ID for ",TRIM(ncvar)
      istatus = 0 
    ENDIF

    ! Looking valtime value
    ncstat = NF_INQ_VARID(ncid,'valtime',vid)
    IF (ncstat .EQ. NF_NOERR) THEN
      ncstat = NF_GET_VAR_INT(ncid,vid,nc_val)
      IF (ncstat .NE. NF_NOERR) THEN
        PRINT *,"Could not find value for Var ID for valtime", vid
        istatus = 0
      ENDIF
    ELSE
      PRINT *,"Could not find Var ID for valtime"
    ENDIF

    ! Looking reftime value
    ncstat = NF_INQ_VARID(ncid,'reftime',vid)
    IF (ncstat .EQ. NF_NOERR) THEN
      ncstat = NF_GET_VAR_INT(ncid,vid,nc_ref)
      IF (ncstat .NE. NF_NOERR) THEN
        PRINT *,"Could not find value for Var ID for reftime", vid
        istatus = 0
      ENDIF
    ELSE
      PRINT *,"Could not find Var ID for reftime"
    ENDIF

    ncstat = NF_CLOSE(ncid)
    RETURN
  END SUBROUTINE ncread_3d
!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ncread_2d(ncfile,nx,ny,ncvar,data2d,nc_val,nc_ref,istatus)

    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN)  :: ncfile
    INTEGER, INTENT(IN)           :: nx,ny
    CHARACTER(LEN=3), INTENT(IN)  :: ncvar
    REAL,INTENT(OUT)              :: data2d(nx,ny)
    INTEGER, INTENT(OUT)          :: istatus

    INTEGER  :: ncid,vid,ncstat
    INTEGER  :: nc_val,nc_ref
    LOGICAL  :: file_exists
    INCLUDE 'netcdf.inc'

    istatus = 1
    INQUIRE(FILE=ncfile,EXIST=file_exists)
    IF (.NOT. file_exists) THEN
      PRINT *, "File not found: ",TRIM(ncfile)
      istatus =0
      RETURN
    ENDIF
    ncstat = NF_OPEN(ncfile,NF_NOWRITE,ncid)
    IF (ncstat .NE. NF_NOERR) THEN
      PRINT *, "Error opening: ",TRIM(ncfile)
      istatus = 0
      RETURN
    ENDIF
    ncstat = NF_INQ_VARID(ncid,ncvar,vid)
    IF (ncstat .EQ. NF_NOERR) THEN
      ncstat = NF_GET_VAR_REAL(ncid,vid,data2d)
      IF (ncstat .NE. NF_NOERR) THEN
        PRINT *, "Problem getting data for ", TRIM(ncvar)
        istatus = 0
      ENDIF
    ELSE
      PRINT *,"Could not find Var ID for ",TRIM(ncvar)
      istatus = 0
    ENDIF

    ! Looking valtime value
    ncstat = NF_INQ_VARID(ncid,'valtime',vid)
    IF (ncstat .EQ. NF_NOERR) THEN
      ncstat = NF_GET_VAR_INT(ncid,vid,nc_val)
      IF (ncstat .NE. NF_NOERR) THEN
        PRINT *,"Could not find value for Var ID for valtime", vid
        istatus = 0
      ENDIF
    ELSE
      PRINT *,"Could not find Var ID for valtime"
    ENDIF

    ! Looking reftime value
    ncstat = NF_INQ_VARID(ncid,'reftime',vid)
    IF (ncstat .EQ. NF_NOERR) THEN
      ncstat = NF_GET_VAR_INT(ncid,vid,nc_ref)
      IF (ncstat .NE. NF_NOERR) THEN
        PRINT *,"Could not find value for Var ID for reftime", vid
        istatus = 0
      ENDIF
    ELSE
      PRINT *,"Could not find Var ID for reftime"
    ENDIF

    ncstat = NF_CLOSE(ncid)
    RETURN
  END SUBROUTINE ncread_2d

!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE calc_val_ref_times(nc_val,nc_ref,time_unit_indicator,time_range,istatus)

    INCLUDE 'lapsparms.cmn' ! laps_cycle_time, model_cycle_time, model_fcst_intvl
    IMPLICIT NONE
    INTEGER, INTENT(OUT)    :: istatus

    INTEGER  :: nc_val, nc_ref, dif
    INTEGER  :: remainder, my_cycle_time
    INTEGER  :: time_unit_indicator, time_range

    istatus = 1
    ! --- Set accum time arg with the use of current time intervals
    dif  = nc_val - nc_ref 
    my_cycle_time = laps_cycle_time
    IF (dif .NE. 0) my_cycle_time = dif

    remainder = MOD(my_cycle_time,3600)
    IF (remainder .EQ. 0) then 
 	! Working with time interval of hours
        time_unit_indicator = 1
        time_range = my_cycle_time / 3600 !Determine no. hrs in laps cycle
        if (dif .EQ. 0) time_range = 0 
    ELSE
 	! Working with time interval of minutes
 	remainder = MOD(my_cycle_time,60)
 	IF (remainder .EQ. 0) then 
           time_unit_indicator = 0
           time_range = my_cycle_time / 60 !Determine no. mins in laps cycle
 	ELSE
           PRINT *,"Problem setting time_unit and time_range"
           istatus = 0
 	ENDIF
    ENDIF

    !PRINT *," ---->> valtime", nc_val, " reftime", nc_ref
    RETURN


  END SUBROUTINE calc_val_ref_times

!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE calc_accum_time(ncfile,ncext,ncvar,nc_val,nc_ref,accum_i4time,time_unit_indicator,startime_accum,&
			     etim_unit,etim_value,istatus)

    INCLUDE 'lapsparms.cmn'     ! laps_cycle_time, model_cycle_time, model_fcst_intvl
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN)  :: ncfile
    CHARACTER(LEN=3), INTENT(IN)  :: ncvar, ncext
    INTEGER, INTENT(OUT)          :: istatus
    INTEGER  :: nc_val, nc_ref, dif
    INTEGER  :: accum_i4time
    CHARACTER(LEN=100) :: str

    INTEGER  :: eyear,emonth,eday,ehour,eminute,esecond
    INTEGER  :: etim_value, etim_unit, remainder, my_cycle_time
    INTEGER  :: time_unit_indicator, startime_accum
    INTEGER  :: ncid,vid,ncstat
    LOGICAL  :: file_exists
    INCLUDE 'netcdf.inc'

    istatus = 1
    IF (ncext .NE. 'fsf') THEN

      INQUIRE(FILE=ncfile,EXIST=file_exists)
      IF (.NOT. file_exists) THEN
        PRINT *, "File not found: ",TRIM(ncfile)
        istatus =0
        RETURN
      ENDIF
      ncstat = NF_OPEN(ncfile,NF_NOWRITE,ncid)
      IF (ncstat .NE. NF_NOERR) THEN
        PRINT *, "Error opening: ",TRIM(ncfile)
        istatus = 0
        RETURN
      ENDIF

      ! Looking for accumulation start time (a9time) in sto_comment, or rto_comment
      ncstat = NF_INQ_VARID(ncid,TRIM(ncvar)//'_comment',vid)
      IF (ncstat .EQ. NF_NOERR) THEN
        ncstat = NF_GET_VAR_TEXT(ncid,vid,str)
        IF (ncstat .NE. NF_NOERR) THEN
          PRINT *,"Could not find value for NC _comment: ", vid
          istatus = 0
        ENDIF
      ELSE
        PRINT *,"Could not find Var ID for ", TRIM(ncvar)//'_comment'
      ENDIF
      ncstat = NF_CLOSE(ncid)

      ! Convert a9time to i4time
      READ(str, '(a9)' ) str
      CALL i4time_fname_lp(str,accum_i4time,istatus)
      CALL cv_i4tim_int_lp(accum_i4time,eyear,emonth,eday,ehour,eminute,esecond,istatus)

      ! --- Set accum time arg with the use of current time intervals
      accum_i4time = accum_i4time - 315619200 
      dif = nc_val - accum_i4time

    ELSE
      !Do not need Var ID for fsf, fua
      dif  = nc_val - nc_ref 
    ENDIF

    my_cycle_time = laps_cycle_time
    IF (dif .NE. 0) my_cycle_time = dif

    remainder = MOD(my_cycle_time,3600)
    IF (remainder .EQ. 0) then 
 	! Working with time interval of hours
 	etim_value = my_cycle_time / 3600
 	etim_unit  = 1 !units of hrs
        time_unit_indicator = 1
        startime_accum = ehour
        IF (dif .EQ. 0) startime_accum = 0 
    ELSE
 	! Working with time interval of minutes
 	remainder = MOD(my_cycle_time,60)
 	IF (remainder .EQ. 0) then 
 	   etim_value = my_cycle_time / 60
 	   etim_unit  = 0 !units of mins
           time_unit_indicator = 0
           startime_accum = eminute
           IF (dif .EQ. 0) startime_accum = 0 
 	ELSE
           PRINT *,"Problem setting accum vars: etim_value and etim_unit"
           istatus = 0
 	ENDIF
    ENDIF

    RETURN

  END SUBROUTINE calc_accum_time
