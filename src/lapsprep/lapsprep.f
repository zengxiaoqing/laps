!dis   
!dis    Open Source License/Disclaimer, Forecast Systems Laboratory
!dis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
!dis    
!dis    This software is distributed under the Open Source Definition,
!dis    which may be found at http://www.opensource.org/osd.html.
!dis    
!dis    In particular, redistribution and use in source and binary forms,
!dis    with or without modification, are permitted provided that the
!dis    following conditions are met:
!dis    
!dis    - Redistributions of source code must retain this notice, this
!dis    list of conditions and the following disclaimer.
!dis    
!dis    - Redistributions in binary form must provide access to this
!dis    notice, this list of conditions and the following disclaimer, and
!dis    the underlying source code.
!dis    
!dis    - All modifications to this software must be clearly documented,
!dis    and are solely the responsibility of the agent making the
!dis    modifications.
!dis    
!dis    - If significant modifications or enhancements are made to this
!dis    software, the FSL Software Policy Manager
!dis    (softwaremgr@fsl.noaa.gov) should be notified.
!dis    
!dis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
!dis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
!dis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
!dis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
!dis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
!dis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
!dis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
!dis   
!dis
 
  PROGRAM lapsprep   
    !
    ! PURPOSE
    ! =======
    ! Prepares LAPS analysis data for ingest by various NWP model pre-processors.
    ! Currently supports MM5V3 (outputs PREGRID v3 format), WRF (outputs 
    ! gribprep format), and RAMS 4.3 (outputs RALPH2 format).   
    !
    ! ARGUMENTS
    ! =========
    !  LAPS_valid_time   - Optional command line argument of form YYJJJHHMM
    !                      specifying time for which to build output
    !                      If not present, the program will use the latest
    !                      available analysis based on LAPS systime.dat
    !
    ! REMARKS
    ! =======
    !  1. You must set the LAPS_DATA_ROOT environment variable before running.
    !  2. Other program controls in lapsprep.nl
    !
    ! HISTORY
    ! =======
    ! 28 Nov 2000 -- Original -- Brent Shaw 
    !    (based on lapsreader program originally developed by Dave Gill of
    !     NCAR to support MM5)

    ! Module declarations

    USE constants
    USE setup
    USE laps_static
    USE lapsprep_mm5
    USE lapsprep_wrf
    USE lapsprep_rams

    ! Variable Declarations

    IMPLICIT NONE

    ! Declarations for use of NetCDF library

    INCLUDE "netcdf.inc" 
    INTEGER :: cdfid , rcode
    INTEGER :: zid
    INTEGER :: z 
    INTEGER , DIMENSION(4) :: start , count
    INTEGER , DIMENSION(2) :: startc, countc
    INTEGER :: vid
    CHARACTER (LEN=132) :: dum 

    ! Arrays for data
    REAL , ALLOCATABLE , DIMENSION (:,:,:) :: u , v , t , rh , ht, &   
                                             lwc,rai,sno,pic,ice, sh, mr, & 
                                             virtual_t, rho
    REAL , ALLOCATABLE , DIMENSION (:,:)   :: slp , psfc, snocov, d2d 
    REAL , ALLOCATABLE , DIMENSION (:)     :: p
    REAL , PARAMETER                       :: tiny = 1.0e-20
    
    ! Miscellaneous local variables
                                        
    INTEGER :: out_loop, loop , var_loop , k, istatus
    LOGICAL :: file_present

    ! Beginning of code
 
    ! Check for command line argument containing LAPS valid time
    ! (in YYJJJHHMM format).  If not present, use the systime.dat
    ! file to get current time.  Note that on HPUX, argument #1
    ! is the executable name and argument #2 is the first actual
    ! argument, unlike other systems.  There is a corrected kludge
    ! for that here.
    CALL GETARG(1,laps_file_time)
    IF (laps_file_time .EQ. 'lapsprep.') THEN 
      ! Must be an HP!
      CALL GETARG(2,laps_file_time)
      IF (laps_file_time .EQ. '         ') THEN
        CALL get_systime(i4time,laps_file_time,istatus)
      ENDIF
    ELSE IF (laps_file_time .EQ. '         ') THEN  
     CALL get_systime(i4time, laps_file_time, istatus)
    ENDIF
    PRINT *, 'LAPS_FILE_TIME = ', laps_file_time
    READ(laps_file_time, '(I2.2,I3.3,I2.2,I2.2)') valid_yyyy, valid_jjj, &
                                                   valid_hh, valid_min
    IF (valid_yyyy.LT.80) THEN
      valid_yyyy = 2000 + valid_yyyy
    ELSE
      valid_yyyy = 1900 + valid_yyyy
    ENDIF
    PRINT '(2A)', 'Running LAPSPREP using A9_time of: ', laps_file_time
  
    ! Get the LAPS_DATA_ROOT from the environment.  

    CALL GETENV('LAPS_DATA_ROOT', laps_data_root)
    PRINT '(2A)', 'LAPS_DATA_ROOT=',laps_data_root

    !  Get the namelist items (from the setup module).

    CALL read_namelist

    ! Get the static information (projection, dimensions,etc.)
 
    PRINT '(A)', 'Getting horizontal grid specs from static file.'
    CALL get_horiz_grid_spec(laps_data_root)

    !  Loop through each of the requested extensions for this date.  Each of the
    !  extensions has a couple of the variables that we want.

    PRINT '(A)', 'Starting Loop for each LAPS file'
    file_loop : DO loop = 1 , num_ext
      !  If this is a microphysical species but not doing 
      !  a hotstart, then cycle over this file.

      IF ( ( (TRIM(ext(loop)).EQ.'lq3').OR.(TRIM(ext(loop)).EQ.'lwc') ).AND. &
          (.NOT.hotstart) ) THEN
        CYCLE file_loop
      ENDIF

      !  Build the input file name.   the input file.

      IF ((TRIM(ext(loop)) .NE. 'lw3' ).AND. &
          (TRIM(ext(loop)) .NE. 'lt1' ).AND. &
          (TRIM(ext(loop)) .NE. 'lh3' )) THEN
        input_laps_file = TRIM(laps_data_root) //'/lapsprd/' // &
            TRIM(ext(loop)) // '/' // laps_file_time // '.' // &
            TRIM(ext(loop))
      ELSE
        IF (balance) THEN
          IF ( ( (TRIM(ext(loop)) .EQ. 'lh3').AND.(adjust_rh) ).OR. &
               (TRIM(ext(loop)) .NE. 'lh3') )THEN
            input_laps_file = TRIM(laps_data_root) //'/lapsprd/balance/' // &
            TRIM(ext(loop)) // '/' // laps_file_time // '.' // &
            TRIM(ext(loop))
          ELSE
            input_laps_file = TRIM(laps_data_root) //'/lapsprd/' // &
              TRIM(ext(loop)) // '/' // laps_file_time // '.' // &
              TRIM(ext(loop))
          ENDIF
        ELSE
          input_laps_file = TRIM(laps_data_root) //'/lapsprd/' // &
              TRIM(ext(loop)) // '/' // laps_file_time // '.' // &
              TRIM(ext(loop)) 
        ENDIF
      ENDIF
      PRINT *, 'Opening: ', input_laps_file

      ! Determine if the file exists

      INQUIRE (FILE=TRIM(input_laps_file), EXIST=file_present)
      IF (.NOT.file_present) THEN
        IF( (ext(loop).EQ.'lt1').OR. &
            (ext(loop).EQ.'lw3').OR. &
            (ext(loop).EQ.'lh3').OR. &
            (ext(loop).EQ.'lsx') ) THEN 
          PRINT '(A)', 'Mandatory file not available!'
          STOP 'not_enough_data'
        ELSE IF ( (ext(loop).EQ.'lq3') .OR. &
               (ext(loop).EQ.'lwc') ) THEN
          PRINT '(A)', 'File not available, cannot do hotstart.'
          hotstart = .false.
          CYCLE file_loop
        ELSE IF (ext(loop).EQ.'l1s') THEN
          PRINT '(A)', 'File not available, cannot do snowcover.'
          CYCLE file_loop
        ELSE
          PRINT '(A)', 'File not available, but not mandatory.'
          CYCLE file_loop
        ENDIF
      ENDIF

      ! Open the netcdf file and get the vertical dimension

      cdfid = NCOPN ( TRIM(input_laps_file) , NCNOWRIT , rcode )

      zid = NCDID ( cdfid , 'z' , rcode )
      CALL NCDINQ ( cdfid , zid , dum , z , rcode )

      IF ( ( ext(loop) .EQ. 'lsx' ) .OR. &
           ( ext(loop) .EQ. 'lm2' ) ) THEN
         z2 = z
      ELSE
         z3 = z
      END IF
      
      IF ( loop .EQ. 1 ) THEN

      ! ALLOCATE space for all of the variables in this data set.  Note
      ! that some ofthe 3d variables are allocated by z+1 instead of just z, as
      ! we are going to put the surface values into these arrays as well.

        ALLOCATE ( u   ( x , y , z3 + 1 ) )
        ALLOCATE ( v   ( x , y , z3 + 1 ) )
        ALLOCATE ( t   ( x , y , z3 + 1 ) )
        ALLOCATE ( rh  ( x , y , z3 + 1 ) )
        ALLOCATE ( ht  ( x , y , z3 + 1 ) )
        ALLOCATE ( slp ( x , y         ) )
        ALLOCATE ( psfc (x , y         ) )
        ALLOCATE ( d2d ( x , y         ) )
        ALLOCATE ( p   (         z3 + 1 ) ) 
        ! The followin variables are not "mandatory"
        ALLOCATE ( lwc ( x , y , z3 ) ) 
        ALLOCATE ( rai ( x , y , z3 ) ) 
        ALLOCATE ( sno ( x , y , z3 ) ) 
        ALLOCATE ( pic ( x , y , z3 ) ) 
        ALLOCATE ( ice ( x , y , z3 ) )
        ALLOCATE ( snocov ( x , y ) ) 
 
        ! The following variables are only used
        ! for converting non-mandatory cloud variables
        ! to mixing ratio values 
        ALLOCATE ( rho ( x , y , z3 ) )
        ALLOCATE ( virtual_t ( x , y , z3 ) )
        ALLOCATE ( sh ( x , y , z ) )
        ALLOCATE ( mr ( x , y , z ) )

        ! Initialize the non-mandatory values

        lwc(:,:,:) = -999.
        rai(:,:,:) = -999.
        sno(:,:,:) = -999.
        pic(:,:,:) = -999.
        snocov(:,:) = -999.
        
      END IF

      IF       ( ext(loop) .EQ. 'lh3' ) THEN

        ! Loop over the number of variables for this data file.

        var_lh3 : DO var_loop = 1 , num_cdf_var(loop)

          ! Get the variable ID.

          vid = NCVID ( cdfid , TRIM(cdf_var_name(var_loop,loop)) , rcode )
          start = (/ 1 , 1 , 1 , 1 /)
          count = (/ x , y , z , 1 /)
          CALL NCVGT ( cdfid , vid , start , count , rh , rcode )

          !  Do this just once for pressure.

          vid = NCVID ( cdfid , 'level' , rcode )

          CALL NCVGT ( cdfid , vid , 1 , z3 , p , rcode )
 
          ! Set the pressure level of the lowest level of our
          ! pressure array as 2001 mb to flag the surface
          p(z3+1) = 2001

        END DO var_lh3 

      ELSE IF (( ext(loop) .EQ. 'lq3' ).AND.(hotstart))THEN

        !  Loop over the number of variables for this data file.
        ! Currently the lq3 file is only used to provide
        ! the specific humidity, which is only used to calculate
        ! virtual temperature so we can calculate density to
        ! scale the microphysical mixing ratios from volume to
        ! dimensionless by mass
        var_lq3 : DO var_loop = 1 , num_cdf_var(loop)

          !  Get the variable ID.

          vid = NCVID ( cdfid , TRIM(cdf_var_name(var_loop,loop)) , rcode )
          start = (/ 1 , 1 , 1 , 1 /)
          count = (/ x , y , z , 1 /)
          CALL NCVGT ( cdfid , vid , start , count , sh , rcode )

        END DO var_lq3
           
      ELSE IF ( ext(loop) .EQ. 'lsx' ) THEN

        !  Loop over the number of variables for this data file.

        var_lsx : DO var_loop = 1 , num_cdf_var(loop)

          !  Get the variable ID.

          vid = NCVID ( cdfid , TRIM(cdf_var_name(var_loop,loop)) , rcode )
          start = (/ 1 , 1 , 1 , 1 /)
          count = (/ x , y , 1 , 1 /)
            
          IF      ( cdf_var_name(var_loop,loop) .EQ. 'u  ' ) THEN
            CALL NCVGT ( cdfid , vid , start , count , u  (1,1,z3+1) , rcode )
          ELSE IF ( cdf_var_name(var_loop,loop) .EQ. 'v  ' ) THEN
            CALL NCVGT ( cdfid , vid , start , count , v  (1,1,z3+1) , rcode )
          ELSE IF ( cdf_var_name(var_loop,loop) .EQ. 't  ' ) THEN
            CALL NCVGT ( cdfid , vid , start , count , t  (1,1,z3+1) , rcode )
          ELSE IF ( cdf_var_name(var_loop,loop) .EQ. 'rh ' ) THEN
            CALL NCVGT ( cdfid , vid , start , count , rh (1,1,z3+1) , rcode )
          ELSE IF ( cdf_var_name(var_loop,loop) .EQ. 'msl' ) THEN
            CALL NCVGT ( cdfid , vid , start , count , slp           , rcode )
          ELSE IF ( cdf_var_name(var_loop,loop) .EQ. 'ps ' ) THEN
            CALL NCVGT ( cdfid , vid , start , count , psfc          , rcode )
          END IF

        END DO var_lsx
      ELSE IF ( ext(loop) .EQ. 'lm2' ) THEN

        var_l1s : DO var_loop = 1 , num_cdf_var(loop)

          !  Get the variable ID.

          vid = NCVID ( cdfid , TRIM(cdf_var_name(var_loop,loop)) , rcode )
          start = (/ 1 , 1 , 1 , 1 /)
          count = (/ x , y , 1 , 1 /)

          IF      ( cdf_var_name(var_loop,loop) .EQ. 'sc ' ) THEN
            CALL NCVGT ( cdfid , vid , start , count , snocov        , rcode )
          END IF

        END DO var_l1s                             
        
      ELSE IF ( ext(loop) .EQ. 'lt1' ) THEN

        !  Loop over the number of variables for this data file.

        var_lt1 : DO var_loop = 1 , num_cdf_var(loop)

          !  Get the variable ID.

          vid = NCVID ( cdfid , TRIM(cdf_var_name(var_loop,loop)) , rcode )
          start = (/ 1 , 1 , 1 , 1 /)
          count = (/ x , y , z , 1 /)
          IF      ( cdf_var_name(var_loop,loop) .EQ. 't3 ' ) THEN
            CALL NCVGT ( cdfid , vid , start , count , t  , rcode )
          ELSE IF ( cdf_var_name(var_loop,loop) .EQ. 'ht ' ) THEN
            CALL NCVGT ( cdfid , vid , start , count , ht , rcode )
          END IF

        END DO var_lt1

      ELSE IF ( ext(loop) .EQ. 'lw3' ) THEN

        !  Loop over the number of variables for this data file.

        var_lw3 : DO var_loop = 1 , num_cdf_var(loop)

          !  Get the variable ID.

          vid = NCVID ( cdfid , TRIM(cdf_var_name(var_loop,loop)) , rcode )
          start = (/ 1 , 1 , 1 , 1 /)
          count = (/ x , y , z , 1 /)
          IF      ( cdf_var_name(var_loop,loop) .EQ. 'u3 ' ) THEN
            CALL NCVGT ( cdfid , vid , start , count , u , rcode )
          ELSE IF ( cdf_var_name(var_loop,loop) .EQ. 'v3 ' ) THEN
            CALL NCVGT ( cdfid , vid , start , count , v , rcode )
          END IF

        END DO var_lw3

      ELSE IF (( ext(loop) .EQ. 'lwc' ).AND.(hotstart)) THEN

        !  Loop over the number of variables for this data file.

        var_lwc1 : DO var_loop = 1 , num_cdf_var(loop)

          !  Get the variable ID.

          vid = NCVID ( cdfid , TRIM(cdf_var_name(var_loop,loop)) , rcode )
          start = (/ 1 , 1 , 1 , 1 /)
          count = (/ x , y , z , 1 /)
          IF      ( cdf_var_name(var_loop,loop) .EQ. 'lwc' ) THEN
            CALL NCVGT ( cdfid , vid , start , count , lwc, rcode )
          ELSE IF ( cdf_var_name(var_loop,loop) .EQ. 'rai' ) THEN
            CALL NCVGT ( cdfid , vid , start , count , rai , rcode )
          ELSE IF ( cdf_var_name(var_loop,loop) .EQ. 'sno' ) THEN
            CALL NCVGT ( cdfid , vid , start , count , sno , rcode ) 
          ELSE IF ( cdf_var_name(var_loop,loop) .EQ. 'ice' ) THEN
            CALL NCVGT ( cdfid , vid , start , count , ice , rcode ) 
          ELSE IF ( cdf_var_name(var_loop,loop) .EQ. 'pic' ) THEN
            CALL NCVGT ( cdfid , vid , start , count , pic , rcode ) 
          END IF

        END DO var_lwc1    

      END IF

    END DO file_loop

    !  Set the lowest level of the geopotential height to topographic height

    ht(:,:,z3+1) = topo 

    IF (hotstart) THEN
      ! If this is a hot start, then we need to convert the microphysical
      ! species from mass per volume to mass per mass (mixing ratio).  This
      ! requires that we compute the air density from virtual temperature
      ! and divide each species by the air density.
       
      ! Convert the specific humidity into mixing ratio.
      mr(:,:,:) = sh(:,:,:) / (1.- sh(:,:,:))
 
      ! Compute virtual temperature from mixing ratio and temperature
      virtual_t(:,:,:)=( 1. + 0.61*mr(:,:,:))*t(:,:,1:z3)
 
      ! Compute density from virtual temperature and gas constant for dry air
      DO k = 1, z3
        rho(:,:,k) = p(k)*100. / (rdry * virtual_t(:,:,k))
      ENDDO
      ! Divide all of species by density     
      lwc(:,:,:) = lwc(:,:,:)/rho(:,:,:)   ! Cloud liquid mixing ratio
      rai(:,:,:) = rai(:,:,:)/rho(:,:,:)   ! Rain mixing ratio
      sno(:,:,:) = sno(:,:,:)/rho(:,:,:)   ! Snow mixing ratio
      ice(:,:,:) = ice(:,:,:)/rho(:,:,:)   ! Ice mixing ratio
      pic(:,:,:) = pic(:,:,:)/rho(:,:,:)   ! Graupel (precipitating ice) mixing rat.
    ENDIF

    ! Loop over the each desired output format

    DO out_loop = 1, num_output
      ! Now it is time to output these arrays.  The arrays are ordered
      !  as (x,y,z).  The origin is the southwest corner at the top of the 
      ! atmosphere for the 3d arrays, where the last layer (z3+1) contains    
      ! the surface information.  This is where you would insert a call
      ! to a custom output routine.

      select_output: SELECT CASE (output_format(out_loop))
        CASE ('mm5 ')
          CALL output_pregrid_format(p, t, ht, u, v, rh, slp, &
                            lwc, rai, sno, ice, pic,snocov)

        CASE ('wrf ')
          CALL output_gribprep_format(p, t, ht, u, v, rh, slp, psfc,&
                             lwc, rai, sno, ice, pic,snocov)
     
        CASE ('rams') 
          CALL output_ralph2_format(p,u,v,t,ht,rh,slp,psfc,snocov)
        CASE ('sfm ')
          PRINT '(A)', 'Support for SFM (RAMS 3b) coming soon...check back later!'

        CASE DEFAULT
          PRINT '(2A)', 'Unrecognized output format: ', output_format
          PRINT '(A)', 'Recognized formats include mm5, rams, wrf, and sfm'

      END SELECT select_output
    ENDDO 
    PRINT '(A)', 'LAPSPREP Complete.'

  END program lapsprep
  


