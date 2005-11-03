      integer maxbgmodels
CWNI-BLS ... 04/19/2005 ... Increased maxbgmodels to 10
CWNI-BLS ... to accomodate UNIDATA netcdf format
CWNI-BLS ... 09/27/2005 ... Increased to 11
CWNI-BLS ... to accomodate WRFARW netCDF data
      parameter (maxbgmodels=11) ! WNI
      real missingflag
      parameter (missingflag=1.0e37)
      real esat(15000:45000), es(15000:45000) 	
      common /estab/ esat,es
