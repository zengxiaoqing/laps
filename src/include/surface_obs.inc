c
c
c.....	LAPS Surface data collection include file.
c......................................................
c.....  This file contains variables for use in the 'obs_driver'
c.....  code that collects the surface data for LAPS.  For porting
c.....  LAPS to different domains, paths to data and some local
c.....  parameters will need to be changed here.  Most grid related
c.....  changes will be done automatically when 'lapsparms.for' is changed.
c
c	include 'lapsparms.for'
c
c.....	Define LAPS grid parameters in terms of 'lapsparms' parameters.
c
c	parameter(ni = NX_L, nj = NY_L, nk = NZ_L)
c
	parameter(badflag = -99.9)
c
c..........................................................................
c..........................................................................
c.....	Changes for local data collection.
c..........................................................................
c..........................................................................
c
c.....	Paths to the data.  Code will have to be changed if there are 
c.....	more than 2 data sources.
c
c	character path_to_METAR*80
c	data path_to_METAR
c     &    /'/public/data/metar/netcdf/'/
c     &    /'/data/fxa/point/metar/netcdf/'/
c
c	character path_to_CDOT*80
c	data path_to_CDOT
c     &    /'/public/data/cdot/netcdf/'/
c     &    /'/data/fxa/point/cdot/'/
c
c.....	Data source (required to set the rest of the path correctly).
c.....
c.....	Set variable to:     1 for FSL /public filename format.
c.....	                     2 for WFO-adv filename format.
c
c	parameter( idata_config = 1 ) 
c
c
c.....	Maximum number of stations for the storage arrays, etc.  Set
c.....	these the same, to the highest number.
c
c	parameter( maxsta = 6000 )	! storage arrays
c	parameter( maxobs = 6000 )	! raw METARs
c
c
c.....	Size of "box" around LAPS grid to collect data for boundary
c.....	analysis.  Around 1 deg lat/lon is good.
c
	parameter( box_size = 1.1 )
c
c
c.....	Time window.  Set these variables for time filtering the data.
c.....	Units are seconds.  For example, if you want data with a time
c.....	stamp from 15 min before to 5 min after the analysis time to
c.....	be included in the LSO file, use 900 and 300 for time_before
c.....	and time_after, respectively.
c
	integer*4 time_before, time_after
	parameter( time_before = 900 )     ! default 15 min.
	parameter( time_after  = 300 )     ! default  5 min.
c
c
c.....	Do we check the METARs for time?  If yes, then off hour LSO
c.....	or LS2 files will not include top of the hour data.  This will
c.....	probably be ok if there are a lot of off hour mesonet data, but
c.....	not so good otherwise.
c
	parameter( ick_METAR_time = 0 )     ! 0=don't check, 1=check
c
c




