      integer maxbgmodels
      parameter (maxbgmodels=13) 
      real missingflag
      parameter (missingflag=1.0e37)
      real esat(15000:45000), es(15000:45000) 	
      common /estab/ esat,es
      character  cbgmodel(maxbgmodels)*10
      data       cbgmodel/'ruc','eta','avn','nfs','fua',' '
     &,' ',' ',' ',' ',' ',' ',' '/
c warning: localization scripting requires the cbgmodel statement above
c All additional comments must be below this line
