!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Program: 	laps2grib
!!
!! Purpose: 	Converts a subset of LAPS analysis grids to a grib2 file
!!
!! Author:	Brent Shaw, Weathernews Inc.
!!		brent.shaw@wni.com
!!
!! History:
!!		7 Dec 2006:	Initial version
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM laps2grib

  !! Requires some modules in LAPS_SRC_ROOT/lib/modules
  USE grib2
  !! Local modules
  USE laps2grib_config
  USE lapsdata
  
  IMPLICIT NONE

  INTEGER                     :: istatus,year,month,day,hour,minute,second
  INTEGER                     :: g2lun
  INTEGER                     :: ncid,varid,ncstat
  INTEGER                     :: k,n
  INTEGER, PARAMETER          :: data_type = 0
  INTEGER, PARAMETER          :: reftime_sig = 0
  INTEGER, PARAMETER          :: process_type = 0
  INTEGER, PARAMETER          :: bg_process_id = 255
  INTEGER, PARAMETER          :: cutoff_hr = 0
  INTEGER, PARAMETER          :: cutoff_min =0 
  INTEGER, PARAMETER          :: time_unit_indicator = 1
  INTEGER, PARAMETER          :: time_range = 0
  INTEGER, PARAMETER          :: levtype_iso = 100
  INTEGER, PARAMETER          :: levscale_iso = 0
  INTEGER, PARAMETER          :: g2miss = 255
  INTEGER, PARAMETER          :: newrec = 1 
  INTEGER, PARAMETER          :: pack_method = 2
  INTEGER                     :: miss_mgmt,inomiss
  
  REAL, ALLOCATABLE           :: lapsdata2d(:,:),lapsdata3d(:,:,:),slab(:,:)
  REAL                        :: r_missing
  CHARACTER(LEN=32),PARAMETER :: vtag = "laps2grib V1.0, 08 Dec 2006, WNI" 
  CHARACTER(LEN=256)          :: laps_data_root,lapsfile
  CHARACTER(LEN=256)          :: g2file, g2file_tmp
  CHARACTER(LEN=512)          :: syscmd

  LOGICAL                     :: lrun_laps2grib

  ! Print banner
  print *, "======================================================"
  print *, "********** ",vtag," **********"
  print *, "======================================================"

  ! Get the LAPS_DATA_ROOT
  CALL GETENV('LAPS_DATA_ROOT',laps_data_root)
  IF ( LEN_TRIM(laps_data_root) .LT. 1) THEN
    PRINT *, "LAPS_DATA_ROOT not set!"
    STOP 'no LAPS_DATA_ROOT'
  ENDIF

  PRINT *, "- LAPS_DATA_ROOT=",TRIM(laps_data_root)
 
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
  CALL get_laps_analtime
  CALL cv_i4tim_int_lp(i4time,year,month,day,hour,minute,second,istatus)
  year = year + 1900

  ! Configure the variable list
  CALL get_data_config(laps_data_root)

  ! Get the rmissing value
  CALL get_r_missing_data(r_missing,istatus)
  PRINT *, "- R_MISSING: ",r_missing
  ! Open the GRIB file
  g2file = TRIM(output_path)//'/'//a9time//'.gr2'
  g2file_tmp = TRIM(output_path)//'/.'//a9time//'.gr2'
  PRINT *, "- Initializing output file: ",TRIM(g2file_tmp)
  CALL init_grib2_file(g2file_tmp,laps_proj,center_id,subcenter_id, &
                       reftime_sig,year,month,day,hour,minute,second, &
                       prod_status,data_type,g2lun,istatus)


  ! Process any 3D isobaric variables
  IF (n_iso .GT. 0) THEN
    ALLOCATE (lapsdata3d(nx,ny,nz))
    ALLOCATE(slab(nx,ny))
    PRINT *, "- Processing 3D Isobaric Fields"
    loop3d: DO n = 1, n_iso

      IF (meta_iso3d(n)%qbal .EQ. 0) THEN
        lapsfile = TRIM(laps_data_root)//'/lapsprd/'//meta_iso3d(n)%ext//'/'//a9time//'.'// &
                         meta_iso3d(n)%ext
      ELSE
        lapsfile = TRIM(laps_data_root)//'/lapsprd/balance/'//meta_iso3d(n)%ext//'/'//&
                     a9time//'.'//meta_iso3d(n)%ext
      ENDIF

      CALL ncread_3d(lapsfile,nx,ny,nz,meta_iso3d(n)%var,lapsdata3d,istatus)
      IF (istatus .EQ. 1) THEN
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
        PRINT *, "Problem getting EXT/VAR",meta_iso3d(n)%ext," ",meta_iso3d(n)%var
      ENDIF
    ENDDO loop3d
    DEALLOCATE(lapsdata3d)
    DEALLOCATE(slab)
  ENDIF

  ! Process any 2D isobaric variables
  IF (n_2d .GT. 0) THEN
    ALLOCATE (lapsdata2d(nx,ny))
    PRINT *, "- Processing 2D Fields"
    loop2d: DO n =1, n_2d
      IF (meta_2d(n)%ext .NE. 'n7g') THEN
        lapsfile = TRIM(laps_data_root)//'/lapsprd/'//meta_2d(n)%ext//'/'//a9time//'.'// &
                         meta_2d(n)%ext

      ELSE
        lapsfile = TRIM(laps_data_root)//'/static/static.nest7grid'
      ENDIF
      CALL ncread_2d(lapsfile,nx,ny,meta_2d(n)%var,lapsdata2d,istatus)
      IF (istatus .EQ. 1) THEN
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
        PRINT *, "Problem getting EXT/VAR",meta_2d(n)%ext," ",meta_2d(n)%var
      ENDIF
    ENDDO loop2d
    DEALLOCATE(lapsdata2d)
  ENDIF
  ! Close the file and rename it to the final output name
  PRINT *, "- Closing File"
  CALL close_grib2_file(g2lun)
  syscmd = 'mv '// TRIM(g2file_tmp) // ' ' // TRIM(g2file)   
  PRINT *, TRIM(syscmd)
  CALL system(syscmd)
  PRINT *, "======================================================"
  PRINT *, "**********      laps2grib completed         **********"
  PRINT *, "======================================================"
END PROGRAM laps2grib

!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ncread_3d(ncfile,nx,ny,nz,ncvar,data3d,istatus)

    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN)  :: ncfile
    INTEGER, INTENT(IN)           :: nx,ny,nz
    CHARACTER(LEN=3), INTENT(IN)  :: ncvar
    REAL,INTENT(OUT)              :: data3d(nx,ny,nz)
    INTEGER, INTENT(OUT)          :: istatus

    INTEGER  :: ncid,vid,ncstat
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
    ncstat = NF_CLOSE(ncid)
    RETURN
  END SUBROUTINE ncread_3d
!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ncread_2d(ncfile,nx,ny,ncvar,data2d,istatus)
                                                                                                  
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN)  :: ncfile
    INTEGER, INTENT(IN)           :: nx,ny
    CHARACTER(LEN=3), INTENT(IN)  :: ncvar
    REAL,INTENT(OUT)              :: data2d(nx,ny)
    INTEGER, INTENT(OUT)          :: istatus
                                                                                                  
    INTEGER  :: ncid,vid,ncstat
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
    ncstat = NF_CLOSE(ncid)
    RETURN
  END SUBROUTINE ncread_2d

