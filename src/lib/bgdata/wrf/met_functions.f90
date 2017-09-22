!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION esat2(t)
    !Computes saturation vapor pressure (Pa) from temperature
    ! NOTE: Computed with respect to liquid!!!
    USE constants
    REAL                       :: esat2
    REAL, INTENT(IN)           :: t
    
    esat2 = 611.21 * EXP ( (17.502 * (t-t0)) / (t-32.18) )
  END FUNCTION esat2 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION mixsat(t,p)
    ! Computes saturation vapor mixing ratio as function of Temp and Press
    USE constants
    IMPLICIT NONE

     REAL,EXTERNAL             :: esat2
     REAL                      :: mixsat
     REAL, INTENT(IN)          :: p
     REAL                      :: satvpr
     REAL, INTENT(IN)          :: t

     satvpr = esat2(t)
     mixsat = (e*satvpr) / (p-satvpr)
   
   END FUNCTION mixsat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION potential_temp(t,p)

  ! Computes potential temperature from temp (K) and pressure (Pa)

    USE constants
    IMPLICIT NONE
    REAL, INTENT(IN)           :: t
    REAL, INTENT(IN)           :: p
    REAL                       :: potential_temp

    potential_temp = t * (p0/p)**kappa

  END FUNCTION potential_temp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION eq_potential_temp(t,p,w,rh)
 
    ! Calculates equivalent potential temperature given temperature, 
    ! pressure, mixing ratio, and relative humidity via
    ! Bolton's equation. (MWR 1980, P 1052, eq. 43)

    IMPLICIT NONE
    
    REAL, INTENT(IN)             :: t  ! temp in K
    REAL, INTENT(IN)             :: p  ! pressure in Pa
    REAL, INTENT(IN)             :: w  ! mixing ratio kg/kg
    REAL, INTENT(IN)             :: rh ! relative humidity (fraction)
    ! REAL, EXTERNAL               :: potential_temp  
    REAL                         :: thtm
    REAL, EXTERNAL               :: tlcl
    REAL                         :: eq_potential_temp

    thtm = t * (100000./p) ** (2. / 7. * ( 1. - (0.28*w)))
    eq_potential_temp = thtm * &
                        EXP ( ( 3.376 / tlcl(t,rh) - 0.00254) * &
                        (w * 1000.0 * (1.0 + 0.81 * w) ) )

  END FUNCTION eq_potential_temp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION tlcl(t,rh)

    ! Computes the temperature of the Lifting Condensation Level
    ! given surface t and RH using Bolton's equation (MWR 1980, p 1048, #22)

    IMPLICIT NONE
    REAL                               :: denom
    REAL, INTENT(IN)                   :: rh
    REAL, INTENT(IN)                   :: t
    REAL                               :: term1
    REAL                               :: term2
    REAL                               :: tlcl

    term1 = 1.0 / (t-55.0)
    term2 = ALOG(rh/1.0)/2840.
    denom = term1 - term2
    tlcl = (1.0/denom) + 55.0

  END FUNCTION tlcl
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION dewpt2(t,rh)
 
  ! Compute dew point from temperature (K) and RH (fraction)
    USE constants
    IMPLICIT NONE
    
    REAL                :: dewpt2
    REAL, INTENT(IN)    :: rh
    REAL, INTENT(IN)    :: t

    dewpt2 = t / (( -rvolv * ALOG(MAX(rh,0.01)) * t) + 1.0 )

  END FUNCTION dewpt2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION relhum(t,mixrat,p)
    ! Computes relative humidity (fraction form)
    IMPLICIT NONE
    REAL, INTENT(IN)                 :: p ! Pressure in Pa
    REAL, INTENT(IN)                 :: mixrat ! vapor mixing ratio
    REAL, EXTERNAL                   :: mixsat ! Saturation vapor mix. ratio
    REAL                             :: relhum  !(fraction)
    REAL, INTENT(IN)                 :: t

    relhum = mixrat/mixsat(t,p)
  END FUNCTION relhum      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION fahren (t)

    ! Converts from Celsius to fahrenheit

    IMPLICIT NONE

    REAL, INTENT(IN)          :: t
    REAL                      :: fahren

    fahren = (1.8 * t) + 32.
  END FUNCTION fahren
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION celsius(tf)

    IMPLICIT NONE
    REAL, INTENT(IN)             :: tf
    REAL                         :: celsius
    ! Converts fahrenheit to celsius

    celsius = (5./9.) * (tf - 32.0)

  END FUNCTION celsius
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION wdir(u,v,xlon,orient,conefact)
    ! Computes wind direction from u/v components.  Converts the direction 
    ! to true direction if conefactor != 0.

    USE constants
    IMPLICIT NONE

    REAL, INTENT(IN)               :: conefact
    REAL                           :: diff
    REAL, INTENT(IN)               :: orient
    REAL, INTENT(IN)               :: u
    REAL, INTENT(IN)               :: v
    REAL                           :: wdir
    REAL, INTENT(IN)               :: xlon

    ! Handle case where u is very small to prevent divide by 0
    
    IF (abs(u) .LT. 0.001) THEN
      IF ( v .le. 0.0) THEN
         wdir = 0.0
      ELSE
         wdir = 180.
      ENDIF

    ! Otherwise, standard trig problem!

    ELSE
      wdir = 270. - (ATAN2(v,u) * rad2deg)
   
      IF  (wdir .gt. 360.) THEN
        wdir = wdir - 360.
      ENDIF
    ENDIF

    ! Change to earth relative

    diff = (orient - xlon) * conefact
    IF (diff .GT. 180.0) diff = diff - 360.
    IF (diff .LT. -180.) diff = diff + 360.

    wdir = wdir - diff
    IF (wdir .gt. 360.0) wdir = wdir - 360.
    IF (wdir .lt. 0.) wdir = wdir + 360.

  END FUNCTION wdir
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION wspd(u,v)

    ! Computes wind velocity from u/v components
    IMPLICIT NONE
    REAL, INTENT(IN)        :: u
    REAL, INTENT(IN)        :: v
    REAL                    :: wspd

    wspd = SQRT ( (u*u) + (v*v) )
  END FUNCTION wspd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION WOBF (T) 

!-----------------------------------------------------------------------
!--   
!--   NAME: WOBUF FUNCTION
!--
!--   THIS FUNCTION CALCULATES THE DIFFERENCE OF THE WET BULB POTENTIAL
!--   TEMPERATURES FOR SATURATED AND DRY AIR GIVEN THE TEMPERATURE.
!--
!--   IT WAS CREATED BY HERMAN WOBUS OF THE NAVY WEATHER RESEARCH
!--   FACILITY FROM DATA IN THE SMITHSONIAN METEOROLOGICAL TABLES.
!--
!-- 
!--   LET WBPTS = WET BULB POTENTIAL TEMPERATURE FOR SATURATED AIR
!--               AT TEMPERATURE T IN CELSIUS
!--
!--   LET WBPTD = WET BULT POTENTIAL TEMPERATURE FOR DRY AIR AT
!--               THE SAME TEMPERATURE.
!--  
!--   THE WOBUS FUNCTION WOBF (IN DEGREES CELSIUS) IS DEFINED BY:
!--
!--               WOBF(T) = WBPTS - WBPTD
!--    
!--   ALTHOUGH WBPTS AND WBPTD ARE FUNCTIONS OF BOTH PRESSURE AND
!--   TEMPERATURE, THEIR DIFFERENCE IS A FUNCTION OF TEMPERATURE ONLY.
!--
!--   THE WOBUS FUNCTION IS USEFUL FOR EVALUATING SEVERAL THERMODYNAMIC
!--   QUANTITIES.
!--
!--   IF T IS AT 1000 MB, THEN T IS POTENTIAL TEMPERATURE PT AND
!--   WBPTS = PT.  THUS,
!--
!--               WOBF(PT) = PT - WBPTD
!--
!--   IF T IS AT THE CONDENSATION LEVEL, THEN T IS THE CONDENSATION
!--   TEMPERATURE TC AND WBPTS IS THE WET BULB POTENTIAL TEMPERATURE
!--   WBPT.  THUS,
!--
!--               WOBF(TC) = WBPT - WBPTD
!--
!--   MANIPULATING THE ABOVE EQUATIONS WE GET,                
!--
!--               WBPT = PT - WOBF(PT) + WOBF(TC)   AND
!--
!--               WBPTS = PT - WOBF(PT) + WOBF(T)
!--
!--   IF T IS EQUIVALENT POTENTIAL TEMPERATURE EPT (IMPLYING THAT
!--   THE AIR AT 1000 MB IS COMPLETELY DRY), THEN 
!--
!--               WBPTS = EPT AND WBPTD = WBPT, THUS,
!--
!--               WOBF(EPT) = EPT - WBPT
!--
!-----------------------------------------------------------------------

      IMPLICIT NONE

      REAL                        :: POL
      REAL,    INTENT(IN)         :: T
      REAL                        :: WOBF
      REAL                        :: X
     
      X = T - 20.0

      IF (X .LE. 0.0) THEN

        POL = 1.0                   + X * (-8.8416605E-03   &
             + X * ( 1.4714143E-04  + X * (-9.6719890E-07   &
             + X * (-3.2607217E-08  + X * (-3.8598073E-10)))))

        WOBF = 15.130 / (POL**4)

      ELSE

        POL = 1.0                   + X * ( 3.6182989E-03   &
             + X * (-1.3603273E-05  + X * ( 4.9618922E-07   &
             + X * (-6.1059365E-09  + X * ( 3.9401551E-11   &
             + X * (-1.2588129E-13  + X * ( 1.6688280E-16)))))))

        WOBF = (29.930 / (POL**4)) + (0.96 * X) - 14.8

      END IF
 
      END FUNCTION WOBF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION heatindex(temp_k, rh_pct)

    ! Computes heat index from a temperature (K) and RH (%).

    USE constants
    IMPLICIT NONE
    REAL, EXTERNAL                  :: celsius
    REAL, EXTERNAL                  :: fahren
    REAL                            :: heatindex
    REAL, INTENT(IN)                :: rh_pct
    REAL                            :: rh_pct_sqr
    REAL                            :: tf
    REAL                            :: tf_sqr
    REAL, INTENT(IN)                :: temp_k

    rh_pct_sqr = rh_pct*rh_pct
    tf = fahren(temp_k - t0)
    tf_sqr = tf*tf

    heatindex =  -42.379 + (2.04901523   * tf)  &
                         + (10.1433312   * rh_pct)  &
                         - (0.22475541   * tf * rh_pct)  &
                         - (6.83783E-03  * tf_sqr)  &
                         - (5.481717E-02 * rh_pct_sqr)   &
                         + (1.22874E-03  * tf_sqr * rh_pct)  &
                         + (8.52E-04     * rh_pct_sqr * tf)  &
                         - (1.99E-06     * tf_sqr * rh_pct_sqr)
 
    heatindex = celsius(heatindex) + t0
  
  END FUNCTION heatindex

  FUNCTION tcvp(p,mr,z,rho,nz)
    USE constants
    IMPLICIT NONE

    REAL                          :: tcvp
    INTEGER, INTENT(IN)           :: nz
    REAL, INTENT(IN)              :: p(nz)
    REAL, INTENT(IN)              :: mr(nz)
    REAL, INTENT(IN)              :: z(nz)
    REAL, INTENT(IN)              :: rho(nz)

    INTEGER   :: k,kbot
    REAL      :: pvapor(nz)
    REAL      :: mrmean, dz

   ! Set top vapor pressure to 0
   pvapor(nz) = 0.

   ! Integrate moisture downward

   DO kbot = nz-1,1,-1
   
     pvapor(kbot) = 0.

     DO k = nz-1,kbot,-1
       
       ! Compute dz and mean Qv for this layer
       dz =z(k+1) - z(k)
       mrmean = (mr(k) + mr(k+1))*0.5

       pvapor(kbot) = pvapor(kbot) + grav*mrmean*rho(k)*dz/(1.+mrmean)
     ENDDO
   ENDDO
   tcvp = pvapor(1)
   RETURN
 END FUNCTION tcvp

    
