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

c       RADAR
        integer*4  MAX_RADAR_FILES
        parameter (MAX_RADAR_FILES = 20000)     

c       Background
        integer*4  MAX_BACKGROUND_FILES
        parameter (MAX_BACKGROUND_FILES = 2000)     

c       Maximum number of LAPS grid levels
        integer*4 MXLVLS
        parameter (MXLVLS=150)

c       Century time cutoff. LAPS will be set to run with filenames having
c       two digits for the year between the actual years of 'iyear_earliest'
c       and 'iyear_earliest+99'. Valid values are from 1901 to 1999. Note that
c       other time constraints may limit the usable time span of LAPS to less
c       than a century.
        integer*4 iyear_earliest
        parameter (iyear_earliest=1950)
