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


MODULE setup

   !  Variables for file names and laps_data_root

   CHARACTER (LEN=256) :: laps_data_root      , &
                          input_laps_file
   CHARACTER (LEN=9)   :: laps_file_time
   INTEGER             :: valid_yyyy, valid_jjj, &
                          valid_hh, valid_min, &
                          i4time

   ! Namelist items

   LOGICAL            :: hotstart,balance,adjust_rh
   CHARACTER (LEN=4)  :: output_format(10)
   INTEGER            :: num_output
   REAL               :: snow_thresh, lwc_scale
   CHARACTER (LEN=256):: output_prefix
  
  !  Output file info.

   INTEGER , PARAMETER :: output = 10
   CHARACTER (LEN=132) :: name

   !  The items below define the number of different LAPS files that will
   !  be read, what their extensions are, and what variables from
   !  each file.  These are hard-coded for now, but we may put this
   !  into lapsprep.nl in the future.  Mandatory extensions must be
   !  listed first

   INTEGER , PARAMETER :: num_ext = 8
   CHARACTER(LEN=3),DIMENSION(num_ext) :: ext = (/ 'lt1' , 'lw3' , &
                                                   'lh3' , 'lsx',  &
                                                   'lsx' , 'lwc',  &
                                                   'lq3' , 'lm2' /)

   CHARACTER(LEN=3),DIMENSION(5,num_ext) :: cdf_var_name = RESHAPE ( &
                           (/ 'ht ' , 't3 ' , 'xxx' , 'xxx' , 'xxx' , & 
                              'u3 ' , 'v3 ' , 'om ' , 'xxx' , 'xxx' , & 
                              'rhl' , 'xxx' , 'xxx' , 'xxx' , 'xxx' , &
                              'u  ' , 'v  ' , 't  ' , 'rh ' , 'tgd' , &
                              'ps ',  'msl' , 'mr ' , 'vv ' , 'xxx' , &
                              'lwc' , 'rai' , 'sno' , 'pic' , 'ice' , &
                              'sh ' , 'xxx' , 'xxx' , 'xxx' , 'xxx' , &
                              'sc ' , 'xxx' , 'xxx' , 'xxx' , 'xxx' /) , &
                                               (/ 5 , num_ext /) )

   INTEGER,DIMENSION(num_ext) :: num_cdf_var = (/2,3,1,5,4,5,1,1/) 

CONTAINS

   SUBROUTINE read_namelist

      IMPLICIT NONE

      INTEGER :: ioerror, nml_unit, c

      NAMELIST /lapsprep_nl/ hotstart      , &
                         balance           , &
                         adjust_rh         , &
                         output_format     , &
                         snow_thresh       , &
                         lwc_scale 

      nml_unit = 77

      output_format(:) = '    '
      ! Set namelist defaults
      balance = .false.
      adjust_rh = .false.
      snow_thresh = 0.5
      lwc_scale = 1.0
 
      ! Open the namelist

      input_laps_file = TRIM(laps_data_root) // '/static/lapsprep.nl'
      OPEN ( FILE = TRIM(input_laps_file)        , &
             UNIT = nml_unit         , &
             STATUS = 'OLD'          , &
             FORM   = 'FORMATTED'    , &
             IOSTAT = ioerror          )

      IF ( ioerror .NE. 0 ) THEN
         PRINT '(3A,I4)','Error opening ', TRIM(input_laps_file), ' Code #',ioerror
         STOP 'error_opening_namelist'
      END IF 

      ! Read and print the namelist

      READ ( nml_unit , NML=lapsprep_nl )
      PRINT '(A)', 'Running LAPSPREP using the following settings:'
      WRITE (6,NML=lapsprep_nl)

      ! Close the namelist

      CLOSE (nml_unit) 
   
      ! Determine number of output formats
      num_output = 0
      find_outputs: DO c = 1 , 10
        IF (output_format(c).EQ.'    ')THEN
          num_output = c - 1
          EXIT find_outputs
        ENDIF
      ENDDO find_outputs
      IF (num_output .LE. 0) THEN
        PRINT '(A)','Must specify at least one output format!'
        STOP 'READ_NAMELIST'
      ENDIF

   END SUBROUTINE read_namelist

END MODULE setup
