  SUBROUTINE saturate_ice_points(sh,t,p,thresh,sh_m,rh_m)

  ! Subroutine to saturate grid boxes with respect to ice 

    IMPLICIT NONE

    ! Inputs:

    REAL, INTENT(IN)    :: sh     ! Specific humidity (kg/kg)   
    REAL, INTENT(IN)    :: t      ! Temperature (K)
    REAL, INTENT(IN)    :: p      ! Pressure (Pa)
    REAL, INTENT(IN)    :: thresh ! Saturation factor    
              ! Set thresh to 1.0 to convert cloud ice up to
              ! ice saturation

    ! Outputs:
   
    REAL, INTENT(OUT)   :: sh_m   ! Adjusted specific humidity
    REAL, INTENT(OUT)   :: rh_m   ! Adjusted RH (%)

    ! Locals

    REAL :: shsat,mr,mrsat,mrmax,tc
    REAL, EXTERNAL :: ssh2,make_rh

    
    ! Set saturation specific humidity for ice for this point         
    tc = t-273.15
    shsat = ssh2(p,tc,tc,0.)*0.001
 
   ! Convert specific humidity to mixing ratio 
    mrsat = shsat/(1.-shsat)
    mr = sh/(1.-sh)
    mrmax = mrsat*thresh
    rh_m = mrmax/mrsat

    ! Convert mrmax to sh_m
    sh_m = mrmax/(1.+mrmax) 
    RETURN
  END SUBROUTINE saturate_ice_points
