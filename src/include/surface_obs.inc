c
c
c.....      LAPS Surface data collection include file.
c......................................................
c.....  This file contains variables for use in the 'obs_driver'
c.....  code that collects the surface data for LAPS.  For porting
c.....  LAPS to different domains, paths to data and some local
c.....  parameters will need to be changed here.  Most grid related
c.....  changes will be done automatically when 'lapsparms.for' is changed.
c
c      include 'lapsparms.for'
c
c.....      Define LAPS grid parameters in terms of 'lapsparms' parameters.
c
c      parameter(ni = NX_L, nj = NY_L, nk = NZ_L)
c
      parameter(badflag = -99.9)
c
c..........................................................................
c..........................................................................
c.....      Changes for local data collection.
c..........................................................................
c..........................................................................
c
c.....      Paths to the data.  Code will have to be changed if there are 
c.....      more than 2 data sources.
c
c      character path_to_METAR*80
c      data path_to_METAR
c     &    /'/public/data/metar/netcdf/'/
c     &    /'/data/fxa/point/metar/netcdf/'/
c
c      character path_to_CDOT*80
c      data path_to_CDOT
c     &    /'/public/data/cdot/netcdf/'/
c     &    /'/data/fxa/point/cdot/'/
c
c.....      Data source (required to set the rest of the path correctly).
c.....
c.....      Set variable to:     1 for FSL /public filename format.
c.....                           2 for WFO-adv filename format.
c
c      parameter( idata_config = 1 ) 
c
c
c.....      Maximum number of stations for the storage arrays, etc.  Set
c.....      these the same, to the highest number.
c
c      parameter( maxsta = 6000 )	! storage arrays
c      parameter( maxobs = 6000 )	! raw METARs
c
c
c.....      Size of "box" around LAPS grid to collect data for boundary
c.....      analysis.  Around 1 deg lat/lon is good.
c
c     parameter( box_size = 1.1 )




