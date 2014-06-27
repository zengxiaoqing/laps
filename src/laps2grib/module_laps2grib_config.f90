!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Module:  laps2grib_config
!!

MODULE laps2grib_config

  ! Modules from the LAPS shared module library
  USE map_utils

  IMPLICIT NONE

  TYPE (proj_info)           :: laps_proj
  INTEGER                    :: nx,ny,nz
  INTEGER                    :: i4time
  INTEGER                    :: num_grids
  CHARACTER(LEN=9)           :: a9time
  REAL, ALLOCATABLE          :: plevels_pa(:)
  CHARACTER(LEN=512)         :: output_path,output_path2
  INTEGER                    :: center_id, subcenter_id
  INTEGER                    :: process_id,prod_status
  LOGICAL                    :: lrun_laps2grib
 
CONTAINS 
  
   SUBROUTINE get_laps_proj(laps_data_root)

     IMPLICIT NONE

     CHARACTER(LEN=*),INTENT(IN) :: laps_data_root
     INTEGER                 :: ncid,dimid,varid,nfstat,proj_code
     CHARACTER(LEN=132)      :: grid_type
     CHARACTER(LEN=256)      :: static_file
     LOGICAL                 :: file_exists
     REAL                    :: dx_m,la1,lo1,latin1,latin2,lov

     INCLUDE 'netcdf.inc'   

     static_file = TRIM(laps_data_root)//'/static/static.nest7grid'
     INQUIRE(FILE=static_file,EXIST=file_exists)
     IF (.NOT.file_exists) THEN
       PRINT *, "-- Static file not found:",TRIM(static_file)
       STOP "no static_file"
     ENDIF

     ! Open the file
     nfstat = NF_OPEN(static_file,NF_NOWRITE,ncid)
     IF (nfstat .NE. NF_NOERR) THEN
       PRINT *, "-- Error returned from NF_OPEN: ",nfstat
       STOP "netcdf file open error"
     ENDIF

     ! Get the data dimensions
     nx = 0
     ny = 0
     nfstat = NF_INQ_DIMID(ncid,'x',dimid)
     nfstat = NF_INQ_DIMLEN(ncid,dimid,nx)
     nfstat = NF_INQ_DIMID(ncid,'y',dimid)
     nfstat = NF_INQ_DIMLEN(ncid,dimid,ny)
     IF ((nx .LT. 1).OR.(ny .LT. 1)) THEN
       PRINT *,"-- Bad nx or ny: ",nx,ny
       STOP "bad dimensions"
     ENDIF

     ! Get the LAPS grid type
     nfstat = NF_INQ_VARID(ncid,'grid_type',varid)
     nfstat = NF_GET_VAR_TEXT(ncid,varid,grid_type)
     IF (grid_type(1:4) .EQ. 'merc') THEN
       proj_code = PROJ_MERC 
     ELSEIF((grid_type(1:4) .EQ. 'seca').OR.&
            (grid_type(1:4) .EQ. 'tang')) THEN 
       proj_code = PROJ_LC
     ELSEIF(grid_type(1:4) .EQ. 'pola') THEN
       proj_code = PROJ_PS
     ELSE
       PRINT *,"-- Unknown LAPS grid_type: ",TRIM(grid_type)
       STOP "bad grid_type"
     ENDIF

     ! Get the grid spacing
     nfstat = NF_INQ_VARID(ncid,'Dx',varid)
     nfstat = NF_GET_VAR_REAL(ncid,varid,dx_m)

     ! Get the truelat1/truelat2/lov
     nfstat = NF_INQ_VARID(ncid,'Latin1',varid)
     nfstat = NF_GET_VAR_REAL(ncid,varid,latin1)
     nfstat = NF_INQ_VARID(ncid,'Latin2',varid)
     nfstat = NF_GET_VAR_REAL(ncid,varid,latin2)
     nfstat = NF_INQ_VARID(ncid,'LoV',varid)
     nfstat = NF_GET_VAR_REAL(ncid,varid,lov)

     ! Get lower left corner lat/lon
     nfstat = NF_INQ_VARID(ncid,'La1',varid)
     nfstat = NF_GET_VAR_REAL(ncid,varid,la1)
     nfstat = NF_INQ_VARID(ncid,'Lo1',varid)
     nfstat = NF_GET_VAR_REAL(ncid,varid,lo1)

     ! Set up the laps_proj structure
     IF (lo1 .GT. 180.) lo1 = lo1 - 360.
     IF (lov .GT. 180.) lov = lov - 360.
     CALL map_set(proj_code,la1,lo1,dx_m,lov,latin1,latin2,nx,ny, &
                  laps_proj)
          
     ! Close the file
     nfstat = NF_CLOSE(ncid)   
     RETURN  
   END SUBROUTINE get_laps_proj
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE get_laps_plevels(laps_data_root)

      IMPLICIT NONE
      
      CHARACTER(LEN=*), INTENT(IN)      :: laps_data_root
      CHARACTER(LEN=256)                :: nl_file
      INTEGER, PARAMETER                :: max_lev = 100
      INTEGER                           :: k
      LOGICAL                           :: file_exists 
      REAL             :: pressures(max_lev)
      NAMELIST /pressures_nl/ pressures

      ! Look for and open the pressures.nl file
      nl_file = TRIM(laps_data_root)//'/static/pressures.nl'
      INQUIRE(FILE=nl_file,EXIST=file_exists)
      IF (.NOT. file_exists) THEN
        PRINT *, "-- No pressures.nl file found: ",TRIM(nl_file)
        STOP "no pressures.nl"
      ENDIF

      pressures(:) = 0.0
      OPEN(FILE=nl_file,UNIT=10,FORM='FORMATTED',STATUS='OLD')
      READ(10,NML=pressures_nl)
      CLOSE(10)

      ! Search through the levels in reverse order, because LAPS
      ! stores the data from the top down and the levels are specified
      ! in pressures_nl from the bottom up

      DO k = max_lev,1,-1
  
        IF (pressures(k) .GT. 0.) THEN

          IF (.NOT. ALLOCATED(plevels_pa)) THEN
            nz = k
            ALLOCATE(plevels_pa(nz))
          ENDIF
          plevels_pa(nz-k+1) = pressures(k)
        ENDIF
      ENDDO
      PRINT *, "-- Number of pressure levels found: ",nz
      PRINT *, plevels_pa 
      RETURN
    END SUBROUTINE get_laps_plevels

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE get_laps_analtime

      IMPLICIT NONE 
      INTEGER   :: istatus

      CALL get_systime(i4time,a9time,istatus)     

      IF (istatus .NE. 1) THEN
        PRINT *, "-- Error getting LAPS analysis time!"
        STOP "get_laps_analtime"
      ENDIF
      PRINT *, "-- LAPS Analysis Time: ", a9time
      RETURN
    END SUBROUTINE get_laps_analtime
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE get_laps_modeltime(modeltime,modeltime_passed_in)
      IMPLICIT NONE 
      CHARACTER(LEN=*), INTENT(IN)   :: modeltime
      INTEGER   :: modeltime_passed_in,istatus

      IF (modeltime_passed_in .EQ. 1) THEN
        print*, "Generating model i4time using ",trim(modeltime)
        CALL i4time_fname_lp(modeltime,i4time,istatus)
        IF (istatus .NE. 1) THEN
          PRINT *, "-- Error converting ",trim(modeltime)," to i4time!"
          STOP "STOP in get_laps_modeltime"
        ENDIF
        a9time = trim(modeltime)
      ELSE
        CALL get_modeltime(i4time,a9time,istatus)
        IF (istatus .NE. 1) THEN
           PRINT *, "-- Error getting LAPS modeltime.dat!"
           STOP "STOP in get_laps_modeltime"
        ENDIF
      ENDIF

      PRINT *, "-- LAPS Model Time: ", a9time
    END SUBROUTINE get_laps_modeltime
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE read_laps2grib_nl(laps_data_root)
  
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN)   :: laps_data_root
      CHARACTER(LEN=256)             :: nl_file
      LOGICAL                        :: file_exists
      NAMELIST /laps2grib_nl/ output_path,output_path2,num_grids,&
                              center_id,subcenter_id,&
                              process_id, prod_status,lrun_laps2grib


      ! Set defaults
      num_grids=1
      output_path(:)  = " "
      output_path2(:) = " "
      center_id = 59 ! NOAA ESRL GSD (aka FSL), for use with templates
      subcenter_id = 1
      process_id = 100
      prod_status = 0
      lrun_laps2grib = .false.

      nl_file = TRIM(laps_data_root)//'/static/laps2grib.nl'
      INQUIRE(FILE=nl_file,EXIST=file_exists)
      IF (file_exists) THEN
        OPEN(FILE=nl_file,UNIT=10,FORM='FORMATTED',STATUS='OLD')
        READ(10,NML=laps2grib_nl)
        CLOSE(10)
      ELSE
        PRINT *, "-- No laps2grib.nl file found: ",TRIM(nl_file)
        PRINT *, "-- Will use default values!"
      ENDIF
      IF (LEN_TRIM(output_path).LT.1) THEN 
        output_path = TRIM(laps_data_root)//'/lapsprd/gr2'
        output_path2 = TRIM(laps_data_root)//'/lapsprd/grb'
      ENDIF
 
      PRINT *, "-- Other configuration: "
      PRINT *, "   Center ID:       ",center_id
      PRINT *, "   Subcenter ID:    ",subcenter_id
      PRINT *, "   Process ID:      ",process_id
      PRINT *, "   Prod Status:     ",prod_status
      PRINT *, "   Output Path:     ",TRIM(output_path)
      RETURN
    END SUBROUTINE read_laps2grib_nl

END MODULE laps2grib_config
