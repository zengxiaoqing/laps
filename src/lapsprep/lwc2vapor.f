  SUBROUTINE lwc2vapor(lwc,sh,t,p,thresh,lwc_m,sh_m,rh_m)

  ! Subroutine to convert cloud water to vapor.

    IMPLICIT NONE

    ! Inputs:

    REAL, INTENT(IN)    :: lwc    ! Cloud water mixing ratio (kg/kg)
    REAL, INTENT(IN)    :: sh     ! Specific humidity (kg/kg)   
    REAL, INTENT(IN)    :: t      ! Temperature (K)
    REAL, INTENT(IN)    :: p      ! Pressure (Pa)
    REAL, INTENT(IN)    :: thresh ! Saturation factor    
              ! Set thresh to 1.0 to convert cloud water up to
              ! vapor saturation.  1.1 will allow 110% RH, and so
              ! forth

    ! Outputs:
   
    REAL, INTENT(OUT)   :: lwc_m  ! Adjusted lwc
    REAL, INTENT(OUT)   :: sh_m   ! Adjusted specific humidity
    REAL, INTENT(OUT)   :: rh_m   ! Adjusted RH (%)

    ! Locals

    REAL :: shsat,shmax
    REAL, EXTERNAL :: ssh,make_rh

    
    ! Set saturation specific humidity for this point         

    shsat = ssh(p,t-273.15)*0.001
    
    shmax = shsat*thresh

    ! Create modified sh (sh_m) by adding cloud liquid

    sh_m = sh + lwc

    ! Zero out the modified cloud water

    lwc_m = 0.

    ! If sh_m exceeds shmax, convert the excess amount 
    ! back to cloud water

    IF (sh_m .GT. shmax) THEN
      lwc_m = sh_m - shmax
      sh_m = shmax
    ENDIF

    !rh_m = (sh_m/shsat)*100.
    rh_m = make_rh(p,t-273.15,sh_m*1000.,-132.0) *100.

    RETURN
  END SUBROUTINE lwc2vapor
