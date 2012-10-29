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


MODULE lapsprep_constants
   REAL , PARAMETER :: pi              = 3.14159265358
   REAL , PARAMETER :: radius_of_earth = 6370. ! MM5/WRF
   REAL , PARAMETER :: radians_per_degree = pi / 180.
   REAL , PARAMETER :: rdry = 287.1  ! Dry air gas constant
   REAL , PARAMETER :: g = 9.81      ! Gravity

   ! Microphysics constants for autoconversion of cloud liquid to
   ! rain and ice to snow, in kg/m**3
   REAL, PARAMETER :: autoconv_lwc2rai = 0.0005
   REAL, PARAMETER :: autoconv_ice2sno = 0.0005
   REAL, PARAMETER :: lwc_min = 0.000001                     
   REAL, PARAMETER :: ice_min = 0.000001                     
   REAL, PARAMETER :: lcp_min = 0.0 ! simpler and more aggressive humidification
END MODULE lapsprep_constants
