cdis    Forecast Systems Laboratory
cdis    NOAA/OAR/ERL/FSL
cdis    325 Broadway
cdis    Boulder, CO     80303
cdis 
cdis    Forecast Research Division
cdis    Local Analysis and Prediction Branch
cdis    LAPS 
cdis 
cdis    This software and its documentation are in the public domain and 
cdis    are furnished "as is."  The United States government, its 
cdis    instrumentalities, officers, employees, and agents make no 
cdis    warranty, express or implied, as to the usefulness of the software 
cdis    and documentation for any purpose.  They assume no responsibility 
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis    
cdis    Permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  All modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making 
cdis    the modifications.  If significant modifications or enhancements 
cdis    are made to this software, the FSL Software Policy Manager  
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 

!       LAPS Parameters
!       Note: Do not include this in the same module as 'lapsparms.cmn'
!             The same variable names are still used in a few cases
!                (i.e. vertical_grid).


!       LAPS Grid Dimensions
!       Notes: The character declaration should equal the length of the string.
!              There should be no trailing blanks in the domain name.
!              The domain name should be in lower case.
        character*9 LAPS_DOMAIN_FILE
        parameter  (LAPS_DOMAIN_FILE = 'nest7grid')

        integer*4 NX_L,NY_L,NZ_L
        parameter (NX_L = 125)
        parameter (NY_L = 105)
        parameter (NZ_L = 21)

!       The following three parameters are also in 'lib/lapsparms.inc'.
        integer*4 NX_L_MAX,NY_L_MAX,NZ_L_MAX
        parameter (NX_L_MAX = 241) ! Used for Local Arrays, must be >= NX_L
        parameter (NY_L_MAX = 241) ! Used for Local Arrays, must be >= NY_L
        parameter (NZ_L_MAX = 41)  ! Used for Local Arrays, must be >= NZ_L

!       This has to do with internal aspects of the wind and temp analyses
        integer*4  IDEN_RATIO_ANAL,IDEN_RATIO_TEMP
        parameter (IDEN_RATIO_ANAL = 6)
        parameter (IDEN_RATIO_TEMP = IDEN_RATIO_ANAL)

        integer*4 I_PERIMETER
        parameter (I_PERIMETER = 10) ! 0

!       parameter VERTICAL_GRID = 'HEIGHT'

        character*8  VERTICAL_GRID
        parameter (VERTICAL_GRID = 'PRESSURE')

        integer*4 NX_M,NY_M
        parameter (NX_M = NX_L)
        parameter (NY_M = NY_L)

        integer*4 laps_cycle_time
        parameter (laps_cycle_time = 3600) ! 1800

        integer*4 i2_missing_data
        parameter (i2_missing_data = -99)
        real*4    r_missing_data
        parameter (r_missing_data = +1e37) ! Also declared in laps_grid_def.h

!       RADAR
        integer*4  MAX_RADARS
        parameter (MAX_RADARS = 3)     ! Also declared in wind/oe/tw_params.h

        real*4 ref_base
        parameter (ref_base = -10.)

        real*4 ref_base_useable
        parameter (ref_base_useable = 0.)    ! >= ref_base

        integer*4  maxstns
        parameter (maxstns = 300) ! For call to read_sfc

        integer*4 N_MESO,N_SAO,N_PIREP
        parameter (N_MESO = maxstns)
        parameter (N_SAO  = maxstns)
        parameter (N_PIREP = 400)

        integer*4 vert_rad_meso
        parameter (vert_rad_meso = 1)  ! Implies vertical spreading of 50mb

        integer*4 vert_rad_sao
        parameter (vert_rad_sao = 1)   ! Implies vertical spreading of 50mb

        integer*4 vert_rad_pirep
        parameter (vert_rad_pirep = 1) ! Implies vertical spreading of 50mb

        integer*4 vert_rad_prof          
        parameter (vert_rad_prof = 1)  ! Implies vertical spreading of 50mb

        integer*4 VAR_LEN, UNITS_LEN, LVL_COORD_LEN
        parameter (VAR_LEN = 3)
        parameter (UNITS_LEN = 10)
        parameter (LVL_COORD_LEN = 4)
 
        integer*4 COMMENT_LEN,DIR_LEN,EXT_LEN
        parameter (COMMENT_LEN = 125)
        parameter (DIR_LEN = 50)
        parameter (EXT_LEN = 31)
 
        integer*4 DOMAIN_NAME_LEN, ASCTIME_LEN
!       DOMAIN_NAME_LEN should be exactly the length of the domain name
        parameter (DOMAIN_NAME_LEN = 9)
        parameter (ASCTIME_LEN = 24)

        integer*4 KDIM_MAX
        parameter (KDIM_MAX = NZ_L_MAX * 6)  ! 6 is max # var in 3D file
!       Following parameter controls the #of sat channels to pss. If this is
!       modified,roce then the Satellite data types parameter below should also be
!       modified.
        integer*4 max_channels 
        parameter (max_channels=5)  !Note this parameter > 0 and <= 5

!       Satellite image processing parameter. Max # of images processed per lvd run.
        integer*4 max_images
        parameter (max_images=2)

!       Radar (vrc) parameters controlling the number of radar data types
!       integer*4 n_radar_types
!        parameter (n_radar_types=1)
!        character*3 c_raddat_types(n_radar_types)
!        data c_raddat_types /'wsi'/ !'wfo'/

!       Parameters for terrain field (generated within static file)
        real*4 silavwt_parm,toptwvl_parm
        parameter (silavwt_parm = 1.)
        parameter (toptwvl_parm = 4.)

c


