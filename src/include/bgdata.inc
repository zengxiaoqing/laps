      integer maxbgmodels
      parameter (maxbgmodels=13) 
      real missingflag
      parameter (missingflag=1.0e37)
      real esat(15000:45000), es(15000:45000) 	
      common /estab/ esat,es
