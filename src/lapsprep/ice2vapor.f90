  SUBROUTINE ice2vapor(ice,sh,t,p,thresh,ice_m,sh_m,rh_m)

  ! Subroutine to convert cloud ice to vapor up to a saturation
  ! threshold (wrt ice)

    IMPLICIT NONE

    ! Inputs:

    REAL, INTENT(IN)    :: ice    ! Cloud ice mixing ratio (kg/kg)
    REAL, INTENT(IN)    :: sh     ! Specific humidity (kg/kg)   
    REAL, INTENT(IN)    :: t      ! Temperature (K)
    REAL, INTENT(IN)    :: p      ! Pressure (Pa)
    REAL, INTENT(IN)    :: thresh ! Saturation factor    
              ! Set thresh to 1.0 to convert cloud ice up to
              ! ice saturation

    ! Outputs:
   
    REAL, INTENT(OUT)   :: ice_m  ! Adjusted lwc
    REAL, INTENT(OUT)   :: sh_m   ! Adjusted specific humidity
    REAL, INTENT(OUT)   :: rh_m   ! Adjusted RH (%)

    ! Locals

    REAL :: shsat,mr,mrsat,mr_m,mrmax,tc
    REAL, EXTERNAL :: ssh2,make_rh

    
    ! Set saturation specific humidity for ice for this point         
    tc = t-273.15
    shsat = ssh2(p,tc,tc,0.)*0.001
   
    ! Convert specific humidity to mixing ratio 
    mrsat = shsat/(1.-shsat)
    mr = sh/(1.-sh)
    mrmax = mrsat*thresh

    ! Create modified mr (mr_m) by adding cloud ice   

    mr_m = mr + ice

    ! Zero out the modified cloud ice  

    ice_m = 0.

    ! If mr_m exceeds mrmax, convert the excess amount 
    ! back to cloud ice   

    IF (mr_m .GT. mrmax) THEN
      ice_m = mr_m - mrmax
      mr_m = mrmax
    ENDIF

    rh_m = (mr_m/mrsat)*100.

    ! Convert mr_m to sh_m
    sh_m = mr_m/(1.+mr_m) 
    RETURN
  END SUBROUTINE ice2vapor
