
      REAL    , PARAMETER :: grav = 9.81  ! acceleration due to gravity (m {s}^-2)

      REAL    , PARAMETER :: r_d          = 287.
      REAL    , PARAMETER :: r_v          = 461.6
      REAL    , PARAMETER :: cp           = 7.*r_d/2.
      REAL    , PARAMETER :: cv           = cp-r_d
      REAL    , PARAMETER :: cliq         = 4190.
      REAL    , PARAMETER :: cice         = 2106.
      REAL    , PARAMETER :: psat         = 610.78
      REAL    , PARAMETER :: rcv          = r_d/cv
      REAL    , PARAMETER :: rcp          = r_d/cp
      REAL    , PARAMETER :: c_2          = cp * rcv
 
      REAL    , PARAMETER :: p1000mb      = 100000.
      REAL    , PARAMETER :: cpovcv       = cp/(cp-r_d)
      REAL    , PARAMETER :: cvovcp       = 1./cpovcv
      REAL    , PARAMETER :: rvovrd       = r_v/r_d
  
      REAL    , PARAMETER :: reradius     = 1./6370.e03

      REAL    , PARAMETER :: asselin      = .025
!     REAL    , PARAMETER :: asselin      = .0
      REAL    , PARAMETER :: cb           = 25.

      REAL    , PARAMETER :: XLV0         = 3.15E6
      REAL    , PARAMETER :: XLV1         = 2370.
      REAL    , PARAMETER :: XLS0         = 2.905E6
      REAL    , PARAMETER :: XLS1         = 259.532

      REAL    , PARAMETER :: XLS          = 2.85E6
      REAL    , PARAMETER :: XLV          = 2.5E6
      REAL    , PARAMETER :: XLF          = 3.50E5

      REAL    , PARAMETER :: rhowater     = 1000.
      REAL    , PARAMETER :: rhosnow      = 100.
      REAL    , PARAMETER :: rhoair0      = 1.28

      REAL    , PARAMETER :: DEGRAD       = 3.1415926/180.
      REAL    , PARAMETER :: DPD          = 360./365.

      REAL    , PARAMETER :: mspkt        = 0.514791
      REAL    , PARAMETER :: ft_per_m     = 1. / (.0254 * 12.)

      REAL    , PARAMETER ::  SVP1=0.6112
      REAL    , PARAMETER ::  SVP2=17.67
      REAL    , PARAMETER ::  SVP3=29.65
      REAL    , PARAMETER ::  SVPT0=273.15
      REAL    , PARAMETER ::  EP_1=R_v/R_d-1.
      REAL    , PARAMETER ::  EP_2=R_d/R_v
      REAL    , PARAMETER ::  KARMAN=0.4
      REAL    , PARAMETER ::  EOMEG=7.2921E-5
      REAL    , PARAMETER ::  STBOLT=5.67051E-8
