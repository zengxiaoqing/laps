      integer maxdprepmodels, maxtimes
      parameter (maxdprepmodels=10)
      parameter (maxtimes=24)
      integer inittimes(maxtimes), fcstlengths(maxtimes)

      character*180 initbgpaths(maxdprepmodels)
     +             ,bndybgpaths(maxdprepmodels)
      integer initbgmodels(maxdprepmodels),bndybgmodels(maxdprepmodels)
      
      namelist /dprep_nl/initbgpaths,bndybgpaths,initbgmodels
     +                  ,bndybgmodels,inittimes, fcstlengths
