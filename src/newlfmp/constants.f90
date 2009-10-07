module constants

   integer, parameter :: imiss=-99999
   integer, parameter :: ismth=0
   integer, parameter :: jsmth=0

   real, parameter :: cm2inch=1.0/2.54
   real, parameter :: cp=1005.7
   real, parameter :: grav=9.81
   real, parameter :: lapse=6.5E-03
   real, parameter :: lv=2.5E+06
   real, parameter :: m2feet=3.281
   real, parameter :: mps2knts=1.944
   real, parameter :: mps2mph=2.237
   real, parameter :: pi=3.141592654
   real, parameter :: p0=100000.0
   real, parameter :: r=287.04         ! 287.058 in Wikipedia
   real, parameter :: rv=461.5
   real, parameter :: t0=273.15
   real, parameter :: xmiss=-99999.9
   real, parameter :: cpog=cp/grav
   real, parameter :: cpor=cp/r
   real, parameter :: deg2rad=pi/180.
   real, parameter :: e=r/rv
   real, parameter :: gor=grav/r
   real, parameter :: kappa=r/cp
   real, parameter :: rad2deg=180./pi
   real, parameter :: rog=r/grav
   real, parameter :: rvolv=rv/lv
   real, parameter :: zero_thresh=1.e-10
   real, parameter :: eradius=6371200.

end module constants
