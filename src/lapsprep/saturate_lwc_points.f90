  SUBROUTINE saturate_lwc_points(sh,t,p,thresh,sh_m,rh_m)

  ! Subroutine to saturate points containing a certain threshold
  ! of cloud liquid

    IMPLICIT NONE

    ! Inputs:

    REAL, INTENT(IN)    :: sh     ! Specific humidity (kg/kg)   
    REAL, INTENT(IN)    :: t      ! Temperature (K)
    REAL, INTENT(IN)    :: p      ! Pressure (Pa)
    REAL, INTENT(IN)    :: thresh ! Saturation factor    
              ! Set thresh to 1.0 to convert cloud water up to
              ! vapor saturation.  1.1 will allow 110% RH, and so
              ! forth

    ! Outputs:
   
    REAL, INTENT(OUT)   :: sh_m   ! Adjusted specific humidity
    REAL, INTENT(OUT)   :: rh_m   ! Adjusted RH (%)

    ! Locals

    REAL :: shsat,mrmax,mr,mr_m,mrsat
    REAL, EXTERNAL :: ssh,make_rh

    
    ! Set saturation specific humidity for this point         

    shsat = ssh(p,t-273.15)*0.001
 
    ! Convert specific humidity to mixing ratio 
    mrsat = shsat/(1.-shsat)
    mr = sh/(1.-sh)
     
    mrmax = mrsat*thresh

    rh_m = thresh*100.

    ! Convert modified mixing ratio back to specific humidity
    sh_m = mrmax/(1.+mrmax)

    RETURN
  END SUBROUTINE saturate_lwc_points 
