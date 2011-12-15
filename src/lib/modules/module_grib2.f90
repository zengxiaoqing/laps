MODULE grib2
  
! Module containing routines to allow output of grib2 data
! Requires the NWS/MDL pk_grib2 library.

  USE map_utils
!  USE constants
  IMPLICIT NONE
  
  PRIVATE 
  INTEGER, PARAMETER                :: kfildo = 6 
  ! Variables to contain GRIB sections
  INTEGER, PARAMETER                :: nidat = 0
  INTEGER, PARAMETER                :: nrdat = 0
  INTEGER                           :: idat ! Normally would have nidat elements
  REAL                              :: rdat ! Normally would have nrdat elements
  INTEGER, PARAMETER                :: ns0 = 16
  INTEGER, PARAMETER                :: ns1 = 21
  INTEGER, PARAMETER                :: ns3 = 96
  INTEGER, PARAMETER                :: ns4 = 60
  INTEGER, PARAMETER                :: ns5 = 49
  INTEGER, PARAMETER                :: ns6 = 6
  INTEGER, PARAMETER                :: ns7 = 8
  INTEGER, PARAMETER                :: ndjer = 30
  INTEGER, PARAMETER                :: l3264b = 32
  INTEGER                           :: idum
  REAL                              :: rdum
  INTEGER, PARAMETER                :: master_table = 1
  INTEGER, PARAMETER                :: local_table = 1
  INTEGER, PARAMETEr                :: grib_edition = 2                                                                                                                                                      
  INTEGER                           :: nd5
  INTEGER                           :: is0 (ns0)
  INTEGER                           :: is1 (ns1)
  INTEGER                           :: is3 (ns3)
  INTEGER                           :: is4 (ns4)
  INTEGER                           :: is5 (ns5)
  INTEGER                           :: is6 (ns6)
  INTEGER                           :: is7 (ns7)
  INTEGER, ALLOCATABLE              :: ib(:,:)
  INTEGER, ALLOCATABLE              :: ipack(:)
  INTEGER                           :: jer(ndjer,2)
  INTEGER                           :: kjer
  INTEGER,PARAMETER                 :: minpk=14
  LOGICAL                           :: big_endian
  INTEGER                           :: ig2status
  LOGICAL                           :: made_sec1
  LOGICAL                           :: made_sec3
 
  PUBLIC open_binary
  PUBLIC init_grib2_file
  PUBLIC close_grib2_file
  PUBLIC close_binary
  PUBLIC write_grib2_template0
  PUBLIC write_grib2_template8
CONTAINS

  SUBROUTINE init_grib2_file(fname,grid,center_id,subcenter_id, &
                             reftime_sig,year,month,day,hour, &
                             minute,second,prod_status,data_type,&
                             funit,istatus)

    ! Opens a new GRIB output file, initializes some of the GRIB
    ! headers
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)  :: fname
    TYPE(proj_info)              :: grid
    INTEGER,INTENT(IN)           :: center_id,subcenter_id
    INTEGER,INTENT(IN)           :: reftime_sig
    INTEGER,INTENT(IN)           :: year,month,day,hour,minute,second
    INTEGER,INTENT(IN)           :: prod_status,data_type
    INTEGER,INTENT(OUT)          :: funit
    INTEGER,INTENT(OUT)          :: istatus

    istatus = 0
    ig2status = 0
    CALL open_binary(fname,funit)
    CALL make_g2_sec1(center_id, subcenter_id, reftime_sig,&
                     year,month,day,hour,minute,second,prod_status,data_type)
    CALL make_g2_sec3(grid)

    IF (ig2status.NE.0) THEN
      PRINT *, "Error Initializing GRIB2 file"
      istatus = 1
    ENDIF
  
  END SUBROUTINE init_grib2_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE close_grib2_file(funit)
  
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: funit
    
    CALL close_binary(funit)
    IF (ALLOCATED(ib)) DEALLOCATE(ib)
    IF (ALLOCATED(ipack)) DEALLOCATE(ipack)
    RETURN
  END SUBROUTINE close_grib2_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE make_g2_sec0(discipline)

    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: discipline  ! Discipline table number
    
    is0(:) = 0
    is0(7) = discipline
    is0(8) = grib_edition
    RETURN
  END SUBROUTINE make_g2_sec0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE make_g2_sec1(center_id, subcenter_id, reftime_sig,&
                          year,month,day,hour,minute,second,prod_status,data_type)

    IMPLICIT NONE
    INTEGER, INTENT(IN)    :: center_id            ! Common code table C-1
    INTEGER, INTENT(IN)    :: subcenter_id         ! 
    INTEGER, INTENT(IN)    :: reftime_sig          ! Signficance of Reference Time, Code Table 1.2, where 0=analysis, 1=fcst
    INTEGER, INTENT(IN)    :: year                 ! Reference time year (4-digits)
    INTEGER, INTENT(IN)    :: month                ! Reference time month
    INTEGER, INTENT(IN)    :: day                  ! Reference time day 
    INTEGER, INTENT(IN)    :: hour                 ! Reference time hour
    INTEGER, INTENT(IN)    :: minute               ! Reference time minute
    INTEGER, INTENT(IN)    :: second               ! Reference time second
    INTEGER, INTENT(IN)    :: prod_status          ! Production status
    INTEGER, INTENT(IN)    :: data_type            ! Code Table 1.4, 0=analysis,1=forecat,2=analysis and forecast
                                                   ! 7=radar observations 

    is1(:) = 255
    is1(5) = 1  ! Section ID
    is1(6) = center_id
    is1(8) = subcenter_id
    is1(10) = master_table
    is1(11) = local_table
    is1(12) = reftime_sig
    is1(13) = year
    is1(15) = month
    is1(16) = day
    is1(17) = hour
    is1(18) = minute
    is1(19) = second
    is1(20) = prod_status
    is1(21) = data_type 
    made_sec1 = .TRUE.
    RETURN
  END SUBROUTINE make_g2_sec1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE make_g2_sec3(grid)

    IMPLICIT NONE
    TYPE(proj_info),INTENT(IN)           :: grid
    REAL                                 :: lat1,lat2,lon1,lon2,lov
    REAL, PARAMETER                      :: llscale = 1.e6
    REAL, PARAMETER                      :: dxscale = 1000.
    REAL, PARAMETER                      :: eradscale = 1
    ig2status = 0 
    ! Initialize is3 

    is3(:) = 0 

    ! Set section number
    is3(5) = 3
    is3(6) = 0  ! Method of specifying grid (0 = uses templates)
    is3(7) = grid%nx * grid%ny
    is3(11) = 0 ! regular grid
    is3(12) = 0 ! regular grid

    ! Starting lat/lon
    lat1 = grid%lat1
    lon1 = grid%lon1
    CALL ij_to_latlon(grid,FLOAT(grid%nx),FLOAT(grid%ny),lat2,lon2)
    lov = grid%stdlon 
    IF (lon1 .LT. 0) lon1 = lon1 + 360.
    IF (lon2 .LT. 0) lon2 = lon2 + 360.
    IF (lov .LT. 0) lov = lov + 360.

    IF (grid%code .EQ. PROJ_MERC) THEN
      is3(13) = 10
      is3(15) = 1  ! Spherical earth, radius provided by user
      is3(16) = 0  ! Scale factor for spherical earth radius
      is3(17) = NINT(earth_radius_m*eradscale)
      is3(31) = grid%nx
      is3(35) = grid%ny
      is3(39) = NINT(lat1 * llscale)
      is3(43) = NINT(lon1 * llscale)
      is3(47) = 56 ! Resolution/component flag
      is3(48) = NINT(grid%truelat1 * llscale)
      is3(52) = NINT(lat2 * llscale)
      is3(56) = NINT(lon2 * llscale)
      is3(60) = 64
      is3(61) = 0
      is3(65) = NINT(grid%dx * dxscale)
      is3(69) = is3(65)
    ELSEIF (grid%code .EQ. PROJ_LC) THEN
      is3(13) = 30
      is3(15) = 1  ! Spherical earth, radius provided by user
      is3(16) = 0  ! Scale factor for spherical earth radius
      is3(17) = NINT(earth_radius_m*eradscale)
      is3(31) = grid%nx
      is3(35) = grid%ny
      is3(39) = NINT(lat1 * llscale)
      is3(43) = NINT(lon1 * llscale)
      is3(47) = 8 ! Resolution/component flag
      is3(48) = NINT(grid%truelat1 * llscale)
      is3(52) = NINT(lov * llscale)
      is3(56) = NINT(grid%dx * dxscale)
      is3(60) = is3(56)
      IF (grid%truelat1 .GT. 0) THEN
        is3(64) = 0
      ELSE
        is3(64) = 128
      ENDIF
      is3(65) = 64
      is3(66) = NINT(grid%truelat1 * llscale)
      is3(70) = NINT(grid%truelat2 * llscale)
      ! I still don't understand the definition of "lat/lon"
      ! of southern pole of project...hope this is correct
      is3(74) = NINT(-90 * llscale)
      is3(78) = NINT(lov * llscale)
                                                                                                                                                        
    ELSEIF (grid%code .EQ. PROJ_PS) THEN
      is3(13) = 20
      is3(15) = 1  ! Spherical earth, radius provided by user
      is3(16) = 0  ! Scale factor for spherical earth radius
      is3(17) = NINT(earth_radius_m*eradscale)
      is3(31) = grid%nx
      is3(35) = grid%ny
      is3(39) = NINT(lat1 * llscale)
      is3(43) = NINT(lon1 * llscale)
      is3(47) = 8 ! Resolution/component flag
      is3(48) = NINT(grid%truelat1 * llscale)
      is3(52) = NINT(lov * llscale)
      is3(56) = NINT(grid%dx * dxscale)
      is3(60) = is3(56)
      IF (grid%truelat1 .GT. 0) THEN
        is3(64) = 0
      ELSE
        is3(64) = 128
      ENDIF
      is3(65) = 64
    ELSE
      PRINT *, "OUTPUT GRIB2:  Projection not supported!"
      ig2status = 1
    ENDIF
   
    ! Allocate space for the packed data
    nd5 = (250 + (grid%nx * grid%ny) + (grid%nx * grid%ny) / 8 + &
         nidat + nrdat )
    IF (ALLOCATED(ipack)) DEALLOCATE(ipack)
    IF (ALLOCATED(ib)) DEALLOCATE(ib)
    ALLOCATE(ipack(nd5))
    ALLOCATE(ib(grid%nx,grid%ny))
    made_sec3 = .TRUE. 
    RETURN
  END SUBROUTINE make_g2_sec3 
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE open_binary(fname,funit)
  
    ! Opens a grib file for writing and returns the funit
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)      :: fname
    INTEGER,INTENT(OUT)                :: funit
    INTEGER, EXTERNAL                  :: c_open_g
    INTEGER                            :: length
    LOGICAL, EXTERNAL                  :: PK_ENDIAN
    length = LEN_TRIM(fname)
    funit = -1
    funit = c_open_g(fname(1:length)//char(0),'w'//char(0))
    big_endian = PK_ENDIAN()
    PRINT *, "OPENED GRIB2 OUTPUT FILE: ", TRIM(fname)
    PRINT *, "  Unit Number: ",funit
    PRINT *, "  BIG_ENDIAN: ",big_endian
    RETURN
  END SUBROUTINE open_binary
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE close_binary(funit)

    IMPLICIT NONE
    INTEGER, INTENT(IN)       :: funit
    INTEGER, EXTERNAL         :: c_close_g
    INTEGER                   :: iretc

    iretc = -1
    iretc = c_close_g(funit)
    print *, 'Grib file closed with iretc = ', iretc
    RETURN
  END SUBROUTINE close_binary
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE write_grib2_template0(funit,discipline,param_category,param_number, &
                                   process_type,bg_process_id, &
                                   af_process_id,cutoff_hr, &
                                   cutoff_min,time_unit_indicator, &
                                   ftime,lev1_type,lev1_scale, lev1_value,&
                                   lev2_type,lev2_scale,lev2_value, &
                                   pack_method,data_scale,miss_mgmt,&
                                   nx,ny,newrec,inomiss,xmissp,xmisss,data_value)
 
    ! This subroutine outputs GRIB2 data to an already opened GRIB2 file, which
    ! was opened by the open_grib2 routine and specified via funit.  
    ! It assumes you have already correctly populated is1 and is3 via
    ! prior calls to make_g2_sec1 and make_g2_sec3.  
    ! This particular routine writes data using Product Definition Template 4.0
    ! and data respresentation template 5.2 (gridpoint complex)
    IMPLICIT NONE

    INTEGER,INTENT(IN)   :: funit     ! File specifier of open GRIB2 file
    INTEGER,INTENT(IN)   :: discipline ! See code table 0.0
    INTEGER,INTENT(IN)   :: param_category  ! Table 4.1
    INTEGER,INTENT(IN)   :: param_number    ! Table 4.2
    INTEGER,INTENT(IN)   :: process_type    ! Table 4.3
    INTEGER,INTENT(IN)   :: bg_process_id   ! Center defined
    INTEGER,INTENT(IN)   :: af_process_id   ! Center defined
    INTEGER,INTENT(IN)   :: cutoff_hr
    INTEGER,INTENT(IN)   :: cutoff_min
    INTEGER,INTENT(IN)   :: time_unit_indicator ! Table 4.4
    INTEGER,INTENT(IN)   :: ftime  ! Forecast time
    INTEGER,INTENT(IN)   :: lev1_type,lev1_scale,lev1_value
    INTEGER,INTENT(IN)   :: lev2_type,lev2_scale,lev2_value
    INTEGER,INTENT(IN)   :: pack_method,data_scale,miss_mgmt
    INTEGER,INTENT(IN)   :: nx,ny,newrec,inomiss
    REAL,INTENT(IN)      :: xmissp,xmisss
    REAL,INTENT(IN)      :: data_value(nx,ny)
    INTEGER              :: dummy1
    INTEGER  :: ibitmap 
    INTEGER  :: imissp, imisss

    imissp = NINT(xmissp)
    imisss = NINT(xmisss)
    ig2status = 0   

    IF (.NOT. made_sec3) THEN
      ig2status = 1
      PRINT *,"Called write_grib2_template0 before section 3 was made!"
      RETURN
    ENDIF
    IF (.NOT. made_sec1) THEN
      ig2status = 1
      PRINT *,"Called write_grib2_template0 before section 1 was made!"
      RETURN
    ENDIF
  
    ! Sections 1 and 3 should already be set up, we ignore section 2,
    ! so finish section 0, 4, 5, and 6 
    ! Set up secion 0 (is0)
    CALL make_g2_sec0(discipline)

    ! Set up section 5 (is5)
    CALL make_g2_sec5(pack_method,data_scale,miss_mgmt)
    IF (ig2status.NE.0) THEN
      PRINT *, "Problem creating Section 5"
      RETURN
    ENDIF

    ! Set up section 6 (bitmap)...currently not using
    is6(:) = 0
    is6(5) = 6
    is6(6) = 255
    ibitmap = 0
  
    ! Set up section 4 ..PDSS
    is4(:) = 0
    is4(5) = 4
    is4(6) = 0
    is4(8) = 0
    is4(10) = param_category
    is4(11) = param_number
    is4(12) = process_type  
    is4(13) = bg_process_id
    is4(14) = af_process_id
    is4(15) = cutoff_hr
    is4(17) = cutoff_min
    is4(18) = time_unit_indicator
    is4(19) = ftime
    is4(23) = lev1_type
    is4(24) = lev1_scale
    is4(25) = lev1_value
    is4(29) = lev2_type
    is4(30) = lev2_scale
    is4(31) = lev2_value

    ! Set up section 7
    is7(:) = 0
    is7(5) = 7
    CALL pk_grib2(kfildo,data_value,dummy1,nx,ny,idat,nidat,rdat,nrdat,&
                  is0,ns0,is1,ns1,is3,ns3,is4,ns4,is5,ns5,is6,ns6, &
                  is7,ns7,ib,ibitmap,ipack,nd5,imissp,xmissp,imisss,xmisss,&
                  newrec,minpk,inomiss,l3264b,jer,ndjer,kjer)

    CALL Error_Check
    IF (.NOT. big_endian) THEN
      CALL c_swap4(nd5*4,ipack)
    ENDIF
    CALL c_write_g(0,is0(9),ipack,funit)

  END SUBROUTINE write_grib2_template0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE make_g2_sec5(pack_method,data_scale,miss_mgmt)

    ! Sets up section 5
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: pack_method,data_scale,miss_mgmt

    is5(:) = 0
    is5(5) = 5

    IF ((pack_method .NE. 0) .AND. (pack_method .NE. 2))THEN

      PRINT *, "Currently, module_grib2 only supports pack_method 0 or 2"
      ig2status = 1
      RETURN
    ENDIF
    is5(10) = pack_method
    is5(18) = data_scale
    IF (pack_method .EQ. 2) THEN
      is5(22) = 1 ! General splitting
      is5(23) = miss_mgmt
    ENDIF
    RETURN
  END SUBROUTINE make_g2_sec5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE write_grib2_template8(funit,discipline,param_category,param_number, &
                                   process_type,bg_process_id, &
                                   af_process_id,cutoff_hr, &
                                   cutoff_min,time_unit_indicator, &
                                   ftime,lev1_type,lev1_scale, lev1_value,&
                                   lev2_type,lev2_scale,lev2_value, &
                                   eyear,emon,eday,ehour,emin,esec,&
                                   ntimes,ntimes_miss,stattype,periodtype,&
                                   etime_unit,etime_value,&
                                   pack_method,data_scale,miss_mgmt,&
                                   nx,ny,newrec,inomiss,xmissp,xmisss,data_value)
                                                                                                                                                                                                                                  
    ! This subroutine outputs GRIB2 data to an already opened GRIB2 file, which
    ! was opened by the open_grib2 routine and specified via funit.
    ! It assumes you have already correctly populated is1 and is3 via
    ! prior calls to make_g2_sec1 and make_g2_sec3.
    ! This particular routine writes data using Product Definition Template 4.8

    IMPLICIT NONE
                                                                                                                                                                                                                                  
    INTEGER,INTENT(IN)   :: funit     ! File specifier of open GRIB2 file
    INTEGER,INTENT(IN)   :: discipline ! See code table 0.0
    INTEGER,INTENT(IN)   :: param_category  ! Table 4.1
    INTEGER,INTENT(IN)   :: param_number    ! Table 4.2
    INTEGER,INTENT(IN)   :: process_type    ! Table 4.3
    INTEGER,INTENT(IN)   :: bg_process_id   ! Center defined
    INTEGER,INTENT(IN)   :: af_process_id   ! Center defined
    INTEGER,INTENT(IN)   :: cutoff_hr
    INTEGER,INTENT(IN)   :: cutoff_min
    INTEGER,INTENT(IN)   :: time_unit_indicator ! Table 4.4
    INTEGER,INTENT(IN)   :: ftime  ! Forecast time
    INTEGER,INTENT(IN)   :: lev1_type,lev1_scale,lev1_value
    INTEGER,INTENT(IN)   :: lev2_type,lev2_scale,lev2_value
    INTEGER,INTENT(IN)   :: eyear,emon,eday,ehour,emin,esec
    INTEGER,INTENT(IN)   :: ntimes,ntimes_miss,stattype,periodtype
    INTEGER,INTENT(IN)   :: etime_unit,etime_value
    INTEGER,INTENT(IN)   :: pack_method,data_scale,miss_mgmt
    INTEGER,INTENT(IN)   :: nx,ny,newrec,inomiss
    REAL, INTENT(IN)     :: xmissp, xmisss
    REAL,INTENT(IN)      :: data_value(nx,ny)
    INTEGER              :: dummy1
    INTEGER  :: ibitmap
    INTEGER  :: imissp, imisss  

    imissp = NINT(xmissp)
    imisss = NINT(xmisss)                                                                                                                                                                                       
    ig2status = 0
                                                                                                                                                                                                                                  
    IF (.NOT. made_sec3) THEN
      ig2status = 1
      PRINT *,"Called write_grib2_template8 before section 3 was made!"
      RETURN
    ENDIF
    IF (.NOT. made_sec1) THEN
      ig2status = 1
      PRINT *,"Called write_grib2_template8 before section 1 was made!"
      RETURN
    ENDIF
                                                                                                                                                                                                                                  
    ! Sections 1 and 3 should already be set up, we ignore section 2,
    ! so finish section 0, 4, 5, and 6
    ! Set up secion 0 (is0)
    CALL make_g2_sec0(discipline)
                                                                                                                                                                                                                                  
    ! Set up section 5 (is5)
    CALL make_g2_sec5(pack_method,data_scale,miss_mgmt)
    IF (ig2status.NE.0) THEN
      PRINT *, "Problem creating Section 5"
      RETURN
    ENDIF
                                                                                                                                                                                                                                  
    ! Set up section 6 (bitmap)...currently not using
    is6(:) = 0
    is6(5) = 6
    is6(6) = 255
    ibitmap = 0

    ! Set up section 4 ..PDSS
    is4(:) = 0
    is4(5) = 4
    is4(6) = 0
    is4(8) = 8
    is4(10) = param_category
    is4(11) = param_number
    is4(12) = process_type
    is4(13) = bg_process_id
    is4(14) = af_process_id
    is4(15) = cutoff_hr
    is4(17) = cutoff_min
    is4(18) = time_unit_indicator
    is4(19) = ftime
    is4(23) = lev1_type
    is4(24) = lev1_scale
    is4(25) = lev1_value
    is4(29) = lev2_type
    is4(30) = lev2_scale
    is4(31) = lev2_value
    ! time of end of overall time interval
    is4(35) = eyear
    is4(37) = emon
    is4(38) = eday
    is4(39) = ehour
    is4(40) = emin
    is4(41) = esec
    ! no. time range specs desc time intvls used to cal statistically-processed field
    is4(42) = ntimes
    is4(43) = ntimes_miss
    is4(47) = stattype
    is4(48) = periodtype
    ! unit of time for time range over which statistical processing is done (Table 4.4)
    is4(49) = etime_unit
    is4(50) = etime_value
    is4(54) = 255
    is4(55:58) = 0

    ! Set up section 7
    is7(:) = 0
    is7(5) = 7
    CALL pk_grib2(kfildo,data_value,dummy1,nx,ny,idat,nidat,rdat,nrdat,&
                  is0,ns0,is1,ns1,is3,ns3,is4,ns4,is5,ns5,is6,ns6, &
                  is7,ns7,ib,ibitmap,ipack,nd5,imissp,xmissp,imisss,xmisss,&
                  newrec,minpk,inomiss,l3264b,jer,ndjer,kjer)
                                                                                                                                                                                                                                  
    IF (.NOT. big_endian) THEN
      CALL c_swap4(nd5*4,ipack)
    ENDIF
    CALL c_write_g(0,is0(9),ipack,funit)
                                                                                                                                                                                                                                  
  END SUBROUTINE write_grib2_template8
 
!#####################################################
   SUBROUTINE Error_Check
 
     INTEGER  :: j

     
     DO j = 1, ndjer

       IF (jer(j,2) .NE. 0) THEN 
         print '("WARNING: GRIB2 WRITE ERROR CODE: ",I4,3x,I2)', jer(j,1),jer(j,2)
       ENDIF
     ENDDO
     RETURN
        


    
   END SUBROUTINE Error_Check
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE grib2
